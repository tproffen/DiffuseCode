module kuplot_math_mod
!
!*****7**************************************************************   
!     Here are all the 'math' routines for KUPLOT.                      
!*****7**************************************************************   
contains
!
SUBROUTINE do_sort (zeile, lp) 
!+                                                                      
!     sort xy file after increasing x values                            
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 1) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
INTEGER                    :: lpara (maxw)
INTEGER  :: ianz, ik, ir, l, i, j 
REAL :: xxx, yyy, dxx, dyy 
REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.ne.1) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
ELSE 
   ik = nint (werte (1) ) 
   IF (ik.lt.iz.and.ik.gt.0) then 
      WRITE (output_io, 1000) ik 
      l = lenc(ik) / 2 + 1 
      ir = lenc(ik) 
!
      loop_1: do
!  10       CONTINUE 
         IF (l.gt.1) then 
            l = l - 1 
            xxx = x (offxy (ik - 1) + l) 
            yyy = y (offxy (ik - 1) + l) 
            dxx = dx (offxy (ik - 1) + l) 
            dyy = dy (offxy (ik - 1) + l) 
         ELSE 
            xxx = x (offxy (ik - 1) + ir) 
            yyy = y (offxy (ik - 1) + ir) 
            dxx = dx (offxy (ik - 1) + ir) 
            dyy = dy (offxy (ik - 1) + ir) 
            x (offxy (ik - 1) + ir) = x (offxy (ik - 1) + 1) 
            y (offxy (ik - 1) + ir) = y (offxy (ik - 1) + 1) 
            dx (offxy (ik - 1) + ir) = dx (offxy (ik - 1) + 1) 
            dy (offxy (ik - 1) + ir) = dy (offxy (ik - 1) + 1) 
            ir = ir - 1 
            IF (ir.eq.1) then 
               x (offxy (ik - 1) + 1) = xxx 
               y (offxy (ik - 1) + 1) = yyy 
               dx (offxy (ik - 1) + 1) = dxx 
               dy (offxy (ik - 1) + 1) = dyy 
!              RETURN 
               exit loop_1
            ENDIF 
         ENDIF 
         i = l 
         j = l + l 
         loop_2: do while(j.le.ir)
!  20       CONTINUE 
!     IF (j.le.ir) then 
            IF (j.lt.ir) then 
               IF (x(offxy(ik - 1) + j) .lt. x(offxy(ik - 1) + j + 1) ) j = j + 1                                  
            ENDIF 
            IF (xxx.lt.x (offxy (ik - 1) + j) ) then 
               x (offxy (ik - 1) + i) = x (offxy (ik - 1) + j) 
               y (offxy (ik - 1) + i) = y (offxy (ik - 1) + j) 
               dx (offxy (ik - 1) + i) = dx (offxy (ik - 1) + j) 
               dy (offxy (ik - 1) + i) = dy (offxy (ik - 1) + j) 
               i = j 
               j = j + j 
            ELSE 
               j = ir + 1 
            ENDIF 
!        GOTO 20 
!     ENDIF 
         end do loop_2
         x (offxy (ik - 1) + i) = xxx 
         y (offxy (ik - 1) + i) = yyy 
         dx (offxy (ik - 1) + i) = dxx 
         dy (offxy (ik - 1) + i) = dyy 
!        GOTO 10 
      enddo loop_1
   ELSE 
      ier_num = - 4 
      ier_typ = ER_APPL 
   ENDIF 
ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Sorting data set ',i3,' ..') 
!
END SUBROUTINE do_sort                        
!
!*****7**************************************************************   
!
SUBROUTINE do_solve (zeile, lp) 
!+                                                                      
!     solve a polynomial for zero points
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
use math_sup
use param_mod
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 5) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lp
!
CHARACTER(LEN=PREC_STRING), dimension(MAXW) :: cpara
INTEGER                   , dimension(MAXW) :: lpara
INTEGER :: ianz
REAL(KIND=PREC_DP)        , dimension(MAXW) :: werte
REAL(KIND=PREC_DP)        , dimension(3   ) :: root
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!
call math_solve_poly(werte(1), werte(2), werte(3), werte(4), ianz, root)
res_para(0) = ianz
res_para(1) = root(1)
res_para(2) = root(2)
res_para(3) = root(3)
!
end SUBROUTINE do_solve
!
!*****7**************************************************************   
!
SUBROUTINE do_four (zeile, lp) 
!+                                                                      
!     Calculates simple discrete Fourier Transform from                 
!     given data set. Range in Fourier space is determined              
!     by sample width in real space.                                    
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
USE precision_mod
USE string_convert_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 5) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
INTEGER :: lpara (maxw) 
INTEGER :: ianz, ik 
REAL(KIND=PREC_DP) :: werte (maxw) 
LOGICAL :: fo (4), lpi 
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.eq.1.or.ianz.eq.2.or.ianz.eq.4.or.ianz.eq.5) then 
   IF (ianz.eq.2.or.ianz.eq.5) then 
      CALL do_cap (cpara (ianz) ) 
      fo(1) = (index(cpara(ianz)(1:lpara(ianz) ), 'R') .ne.0)
      fo(2) = (index(cpara(ianz)(1:lpara(ianz) ), 'I') .ne.0)
      fo(3) = (index(cpara(ianz)(1:lpara(ianz) ), 'A') .ne.0)
      fo(4) = (index(cpara(ianz)(1:lpara(ianz) ), 'P') .ne.0)
      lpi = (index(cpara(ianz)(1:lpara(ianz) ), 'H') .ne.0) 
      ianz = ianz - 1 
!                                                                       
      IF (.not. (fo (1) .or.fo (2) .or.fo (3) .or.fo (4) ) ) then 
         fo (1) = .true. 
      ENDIF 
   ELSE 
      fo (1) = .true. 
      fo (2) = .false. 
      fo (3) = .false. 
      fo (4) = .false. 
      lpi = .false. 
   ENDIF 
!                                                                       
   CALL ber_params (ianz, cpara, lpara, werte, maxw) 
   ik = nint (werte (1) ) 
   IF (ik.lt.iz.and.ik.gt.0) then 
      IF (lni (ik) ) then 
         IF (ianz.le.2) then 
            CALL do_four_z (ik, fo, lpi) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         CALL do_four_y (ik, fo, lpi, ianz, werte, maxw) 
      ENDIF 
   ELSE 
      ier_num = - 4 
      ier_typ = ER_APPL 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
END SUBROUTINE do_four                        
!
!*****7**************************************************************** 
!
SUBROUTINE do_four_y (ik, fo, lpi, ianz, werte, maxw) 
!                                                                       
!     This routine perfoms the Fourier Transform for 2D data            
!                                                                       
USE errlist_mod 
USE prompt_mod 
USE wink_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
integer                            , intent(in) :: ik
logical, dimension(4)              , intent(in) :: fo
logical                            , intent(in) :: lpi
integer                            , intent(in) :: ianz
integer                            , intent(in) :: MAXW
REAL(KIND=PREC_DP), dimension(MAXW), intent(in) :: werte (maxw) 
!
CHARACTER(len=10) :: oname (4) 
INTEGER :: i, j, ipkt, inew, k, kk 
INTEGER :: maxpkt 
INTEGER :: ift (0:4) 
REAL(kind=PREC_DP) :: dxx, h, hmin, hmax, dh, fc, f (4) 
!     LOGICAL fo (4), lpi 
!                                                                       
DATA oname / 'real.xy', 'imag.xy', 'ampl.xy', 'phase.xy' / 
!                                                                       
!------ Enough space for new data sets                                  
!                                                                       
inew = 0 
ift (0) = ik 
DO i = 1, 4 
   IF (fo (i) ) then 
      inew = inew + 1 
      ift (i) = iz - 1 + inew 
   ENDIF 
ENDDO 
!                                                                       
IF (inew.eq.0) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
IF ( (iz + inew - 1) .gt.maxkurvtot) then 
   ier_num = - 1 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ Figure out range                                                
!                                                                       
dxx = x (offxy (ik - 1) + 2) - x (offxy (ik - 1) + 1) 
IF (ianz.eq.4) then 
   hmin = werte (2) 
   hmax = werte (3) 
   dh = werte (4) 
ELSE 
   hmin = - lenc(ik) / dxx 
   hmax = lenc(ik) / dxx 
   dh = 1.0 / (lenc(ik) * dxx) 
ENDIF 
ipkt = nint ( (hmax - hmin) / dh) + 1 
!                                                                       
maxpkt = maxarray - offxy (iz - 1) 
IF (inew * ipkt.gt.MAXARRAY) then 
   ier_num = - 6 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ Setting properties of newly created data sets                   
!                                                                       
DO i = 1, 4 
   IF (fo (i) ) then 
      fname (ift (i) ) = oname (i) 
      fform (ift (i) ) = 'XY' 
      lni (ift (i) ) = .false. 
      lenc(ift (i) ) = ipkt 
      offxy (ift (i) ) = offxy (ift (i) - 1) + lenc(ift (i) ) 
      offz (ift (i) ) = offz (ift (i) - 1) 
      iz = iz + 1 
   ENDIF 
ENDDO 
!                                                                       
!------ Calculating the Fourier Transform                               
!                                                                       
IF (lpi) then 
   fc = zpi 
ELSE 
   fc = 1. 
ENDIF 
!                                                                       
WRITE (output_io, 1000) ik, lpi 
!                                                                       
kk = 0 
DO k = 1, ipkt 
   h = hmin + dh * (k - 1) 
   f (1) = 0.0 
   f (2) = 0.0 
   DO j = 1, lenc(ik) 
      f (1) = f(1) + y(offxy(ik - 1) + j) * cos(fc * h * x(offxy(ik - 1) + j) )
      f (2) = f(2) + y(offxy(ik - 1) + j) * sin(fc * h * x(offxy(ik - 1) + j) )
   ENDDO 
   f (1) = dxx * f (1) 
   f (2) = dxx * f (2) 
   f (3) = sqrt (f (1) **2 + f (2) **2) 
   f (4) = atan2(f (2) , f (1) ) / rad 
   kk = kk + 1 
!                                                                       
   DO i = 1, 4 
      IF (fo (i) ) then 
         x (offxy (ift (i) - 1) + kk) = h 
         y (offxy (ift (i) - 1) + kk) = f (i) 
         dx (offxy (ift (i) - 1) + kk) = 0.0 
         dy (offxy (ift (i) - 1) + kk) = 0.0 
      ENDIF 
   ENDDO 
ENDDO 
!                                                                       
CALL get_extrema 
!                                                                       
 1000 FORMAT     (' ------ > Calulating Fourier Transform of data set ',&
     &                              I2,' (2pi=',L1,') ...')             
END SUBROUTINE do_four_y                      
!
!*****7**************************************************************** 
!
SUBROUTINE do_four_z (ik, fo, lpi) 
!                                                                       
!     This routine perfoms the Fourier Transform for 3D data            
!                                                                       
USE errlist_mod 
USE prompt_mod 
USE wink_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
integer                            , intent(in) :: ik
logical, dimension(4)              , intent(in) :: fo
logical                            , intent(in) :: lpi
!
CHARACTER(len=12) :: oname (4) 
INTEGER :: i, jx, jy, ipkt, jpkt, inew, kkx, kky 
INTEGER :: kx, ky, maxpkt, maxzz 
INTEGER :: ift (0:4) 
REAL(kind=PREC_DP) :: dxx, dyy, hx, hy, fc, f (4) 
!                                                                       
DATA oname / 'real.nipl', 'imag.nipl', 'ampl.nipl', 'phase.nipl' / 
!                                                                       
!------ Enough space for new data sets                                  
!                                                                       
inew = 0 
ift (0) = ik 
DO i = 1, 4 
   IF (fo (i) ) then 
      inew = inew + 1 
      ift (i) = iz - 1 + inew 
   ENDIF 
ENDDO 
!                                                                       
IF (inew.eq.0) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
IF ( (iz + inew - 1) .gt.maxkurvtot) then 
   ier_num = - 1 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
maxpkt = maxarray - offxy (iz - 1) 
maxzz = maxarray - offz (iz - 1) 
!                                                                       
ipkt = 2 * (nx (ik) / 2) + 1 
jpkt = 2 * (ny (ik) / 2) + 1 
IF (nx (ik) .le.0.or.ny (ik) .le.0.or.inew * nx (ik) * ny (ik)    &
      .gt.maxzz.or.inew * max (nx (ik), ny (ik) ) .gt.maxpkt) then      
   ier_num = - 6 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ Setting properties of newly created data sets                   
!                                                                       
DO i = 1, 4 
   IF (fo (i) ) then 
      fname (ift (i) ) = oname (i) 
      fform (ift (i) ) = 'NI' 
      lni (ift (i) ) = .true. 
      nx (ift (i) ) = ipkt 
      ny (ift (i) ) = jpkt 
      lenc(ift (i) ) = max (ipkt, jpkt) 
      offxy (ift (i) ) = offxy (ift (i) - 1) + lenc(ift (i) ) 
      offz (ift (i) ) = offz (ift (i) - 1) + nx (ift (i) ) * ny (ift (i) )
      iz = iz + 1 
   ENDIF 
ENDDO 
!                                                                       
!------ Calculating the Fourier Transform                               
!                                                                       
WRITE (output_io, 1000) ik 
!                                                                       
IF (lpi) then 
   fc = REAL(zpi) 
ELSE 
   fc = 1. 
ENDIF 
!                                                                       
kkx = 0 
dxx = x (offxy (ik - 1) + 2) - x (offxy (ik - 1) + 1) 
dyy = y (offxy (ik - 1) + 2) - y (offxy (ik - 1) + 1) 
DO kx = - nx (ik) / 2, nx (ik) / 2 
   kkx = kkx + 1 
   kky = 0 
   hx = kx / (nx (ik) * dxx) 
   DO ky = - ny (ik) / 2, ny (ik) / 2 
      kky = kky + 1 
      hy = ky / (ny (ik) * dyy) 
      f (1) = 0.0 
      f (2) = 0.0 
      DO jx = 1, nx (ik) 
         DO jy = 1, ny (ik) 
            f(1) = f(1) + z(offz(ik - 1) + (jx - 1) * ny(ik) + jy)       &
                   * cos(fc * (hx * x (offxy(ik - 1) + jx) + hy * y(offxy(ik - 1) + jy) ) )
            f(2) = f(2) + z(offz(ik - 1) + (jx - 1) * ny(ik) + jy)       &
                   * sin(fc * (hx * x (offxy(ik - 1) + jx) + hy * y(offxy(ik - 1) + jy) ) )
         ENDDO 
      ENDDO 
      f (1) = dxx * dyy * f (1) 
      f (2) = dxx * dyy * f (2) 
      f (3) = sqrt (f (1) **2 + f (2) **2) 
      f (4) = atan2(f (2) , f (1) ) / REAL(rad) 
!                                                                       
      DO i = 1, 4 
         IF (fo (i) ) then 
            z (offz (ift (i) - 1) + (kkx - 1) * ny (ift (i) ) + kky)       &
            = f (i)                                                        
            x (offxy (ift (i) - 1) + kkx) = hx 
            y (offxy (ift (i) - 1) + kky) = hy 
            dx (offxy (ift (i) - 1) + kkx) = 0.0 
            dy (offxy (ift (i) - 1) + kky) = 0.0 
         ENDIF 
      ENDDO 
   ENDDO 
ENDDO 
!                                                                       
CALL get_extrema 
!                                                                       
 1000 FORMAT     (' ------ > Calulating Fourier Transform of data set ',&
     &                              I2,' ...')                          
!
END SUBROUTINE do_four_z                      
!
!!*****7**************************************************************   
!
SUBROUTINE kuplot_do_fft(zeile, lp) 
!-
!  Calculate a FFT of data set(s)
!+
!
USE kuplot_mod
use kuplot_global
!
USE errlist_mod
USE get_params_mod
use lib_math_mod
USE map_1dtofield
USE precision_mod
use string_convert_mod
USE take_param_mod
use wink_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER         , INTENT(INOUT) :: lp
!
!
INTEGER, PARAMETER :: MAXW = 4     ! Max parameters == NOPTIONAL
INTEGER, PARAMETER :: MAXU = 3     ! Max parameters at O_PAD
CHARACTER(LEN=PREC_STRING)                  :: string
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER                   , DIMENSION(MAXW) :: lpara
INTEGER :: ianz
INTEGER :: jjanz
INTEGER :: isdim                   ! Dimension of data set (1, 2, 3)
!
INTEGER, DIMENSION(2) :: idata     ! Data set number to transform
INTEGER, DIMENSION(2) :: odata     ! Data set number with results
integer, dimension(3) :: ipad      ! Pad data to this dimension with zeros outside
real(kind=PREC_DP) :: xscale
real(kind=PREC_DP) :: yscale
REAL(KIND=PREC_DP) , DIMENSION(MAXU) :: uwerte   ! Calculated values
!
INTEGER, PARAMETER :: NOPTIONAL = 4
INTEGER, PARAMETER :: O_REAL      = 1
INTEGER, PARAMETER :: O_IMAG      = 2
INTEGER, PARAMETER :: O_PAD       = 3
INTEGER, PARAMETER :: O_STYLE     = 4
!INTEGER, PARAMETER :: O_COLDY     = 5
!INTEGER, PARAMETER :: O_LAYER     = 6
!INTEGER, PARAMETER :: O_SEPARATOR = 7
!INTEGER, PARAMETER :: O_DECIMAL   = 8
CHARACTER(LEN=          5), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 2 ! Number of values to calculate 
!
DATA oname  / 'real ', 'imag ' , 'pad'   , 'style' /
DATA loname /  4     ,  4      ,  3      ,  5      /
opara  =  (/ '1.00000', '1.00000', '[0,0,0]', 'Q      ' /)
lopara =  (/  6,         6       ,  7       ,  1       /)
owerte =  (/  1.0,       1.0     ,  0.0     ,  0.0     /)
idata  = 0
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp)
IF(ier_num /= 0) RETURN
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_num /= 0) RETURN
!
string = opara(O_PAD)
lp     = lopara(O_PAD)
jjanz  = 3
call get_optional_multi(MAXU, string, lp, uwerte, jjanz)
if(ier_num/=0) return
if(minval(uwerte)<0.0D0) then
   ier_num = -77
   ier_typ = ER_APPL
   ier_msg(1) = 'Check values at pad:[ , , ] parameter'
   return
endif
ipad = nint(uwerte)
!
call do_cap(opara(O_STYLE))
!
!
isdim = 1
IF(lpresent(O_REAL) .AND. lpresent(O_IMAG)) THEN        ! User specified REAL and IMAG
   idata(1) = NINT(owerte(O_REAL))
   idata(2) = NINT(owerte(O_IMAG))
   IF((lh5(idata(1)).EQV.lh5(idata(2)))  .AND.   &        ! Both DATA equivalent ?
      (lni(idata(1)).EQV.lni(idata(2)))        ) THEN     ! Both DATA equivalent ?
      IF(lh5(idata(1))) THEN                            ! Both are 3D
         isdim = ku_ndims(idata(1))                     ! Use actual dimensions
      ELSEIF(lni(idata(1))) THEN                        ! Both are 2D
         isdim = 2
         IF(.NOT. (nx(idata(1))==nx(idata(2)) .AND.    &
                   ny(idata(1))==ny(idata(2))  )) THEN  ! Different length
            ier_num = -74
            ier_typ = ER_APPL
            RETURN
         ENDIF
      ELSE                                              ! Both are 1D
         isdim = 1    
         IF(.NOT. (lenc(idata(1))==lenc(idata(2)))) THEN   ! Different length
            ier_num = -74
            ier_typ = ER_APPL
            RETURN
         ENDIF
      ENDIF
   ELSE
      ier_num = -73                                     ! Data sets differ in dimension
      ier_typ = ER_APPL
      RETURN
   ENDIF
ELSEIF(lpresent(O_REAL) .AND. .NOT.lpresent(O_IMAG)) THEN  ! Only REAL part present
   idata(1) = NINT(owerte(O_REAL))
   idata(2) = 0
   IF(lh5(idata(1))) THEN                                  ! 3D data set
      isdim = ku_ndims(idata(1))                           ! Use actual dimensions
   ELSEIF(lni(idata(1))) THEN                              ! 2D data set
      isdim = 2
   ELSE
      isdim = 1
   ENDIF
ELSEIF(.NOT.lpresent(O_REAL) .AND. lpresent(O_IMAG)) THEN  ! Only IMAG part present
   idata(1) = 0
   idata(2) = NINT(owerte(O_REAL))
   IF(lh5(idata(2))) THEN                                  ! 3D data set
      isdim = ku_ndims(idata(1))                           ! Use actual dimensions
   ELSEIF(lni(idata(2))) THEN                              ! 2D data set
      isdim = 2
   ELSE
      isdim = 1
   ENDIF
ELSE
   ier_num = -73
   ier_typ = ER_APPL
   RETURN
ENDIF
!
xscale = ZPI
yscale = 1.0D0/sqrt(ZPI)
if(opara(O_STYLE) == 'Q') then
   xscale = ZPI
   yscale = 1.0D0/(sqrt(ZPI))**isdim
elseif(opara(O_STYLE) == 'H') then
   xscale = 1.0D0
   yscale = 1.0D0
endif
!
lh5(0) = .true.
IF(isdim==3) THEN                                          ! 2D FFT
!  CALL kuplot_do_fft_3D(idata, odata, xscale, yscale, ipad)
   odata(1) = iz
   odata(2) = iz+1
   CALL fft_3D_global(idata, odata, xscale, yscale, ipad)
   call data3nipl(odata(1), 'Fourier Real', .true.)
   call data3nipl(odata(2), 'Fourier Imag', .true.)
ELSEIF(isdim==2) THEN                                      ! 2D FFT
   if(lh5(idata(1)) .and. lh5(idata(2))) then
!     CALL kuplot_do_fft_2D_global_n(idata, odata, xscale, yscale, ipad)
!
      odata(1) = iz
      odata(2) = iz+1
      call fft_2D_global(idata, odata, xscale, yscale, ipad)    ! Perform FFT on global data
      call data2nipl(odata(1), 'Fourier Real', .true.)
      call data2nipl(odata(2), 'Fourier Imag', .true.)
   else
      CALL kuplot_do_fft_2D(idata, odata, xscale, yscale)
   endif
ELSEIF(isdim==1) THEN                                      ! 1D FFT
   if(lh5(idata(1)) .and. lh5(idata(2))) then
!     CALL kuplot_do_fft_1D_global(idata, odata, xscale, yscale, ipad)
      odata(1) = iz
      odata(2) = iz+1
      call fft_1D_global(idata, odata, xscale, yscale, ipad)    ! Perform FFT on global data
      call data2line(odata(1), 'Fourier Real', .true.)
      call data2line(odata(2), 'Fourier Imag', .true.)
   else
      CALL kuplot_do_fft_1D(idata, odata, xscale, yscale)
   endif
ENDIF
lh5(0) = .false.
!
END SUBROUTINE kuplot_do_fft
!
!*******************************************************************************
!
SUBROUTINE kuplot_do_fft_1D(idata, odata, xscale, yscale)
!-
!  Calculate 1D FFT of data set idata, use the local data structure
!+
!
USE kuplot_mod
!use kuplot_global
use kuplot_show_mod
!
use errlist_mod
USE map_1dtofield
USE wink_mod
use lib_f90_fftw3
!
use iso_c_binding
!
IMPLICIT NONE
!
INTEGER, DIMENSION(2), INTENT(IN)  :: idata     ! Data set number to transform
INTEGER, DIMENSION(2), INTENT(OUT) :: odata     ! Data set number to transform
real(kind=PREC_DP)   , intent(in)  :: xscale
real(kind=PREC_DP)   , intent(in)  :: yscale
!
INTEGER :: i
INTEGER :: kdat
INTEGER :: length             ! Data set length
INTEGER, DIMENSION(3) :: num  ! DATA set dimensions
INTEGER, DIMENSION(3) :: dsort
COMPLEX(KIND=KIND(0.0D0)) , DIMENSION(:), ALLOCATABLE  :: k_data   ! The Kuplot in data set)
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:), ALLOCATABLE  ::  in_pattern  ! The Curve to be FFT'd
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:), ALLOCATABLE  :: out_pattern  ! The Result   of FFT
REAL :: xrange
REAL :: xstep
type(c_ptr) :: plan    ! FFWT3 plan
!
!write(*,*) ' IDATA(1)',     (idata(1))
!write(*,*) ' LENC(1) ', lenc(idata(1)), ' :', lenc(lbound(lenc,1):4)
!write(*,*) ' OFFXY   ', offxy(0:2)
!write(*,*) ' IZ      ', iz
!
IF(idata(1)>0) THEN
   length = lenc(idata(1))                      ! User provided real part
   kdat   = 1
ELSE
   length = lenc(idata(2))                      ! User provided imag part
   kdat   = 2
ENDIF
num    = 1
num(1) = length
dsort(1) = 1
dsort(2) = 2
dsort(3) = 3
ALLOCATE(k_data (length))
ALLOCATE( in_pattern(length))
ALLOCATE(out_pattern(length))
IF(idata(1)>0 .and. idata(2)>0) THEN            ! Got real and imag part
   DO i=1,length
      k_data(i) = CMPLX(y(offxy(idata(1)-1)+i), y(offxy(idata(2)-1)+i))
   ENDDO
ELSEIF(idata(1)>0) THEN                         ! Got real only
   DO i=1,length
      k_data(i) = CMPLX(y(offxy(idata(1)-1)+i), 0.0D0)
   ENDDO
ELSEIF(idata(2)>0) THEN                         ! Got imag only
   DO i=1,length
      k_data(i) = CMPLX(0.0D0, y(offxy(idata(2)-1)+i))
   ENDDO
ENDIF
!
CALL maptofftfd(num, dsort, k_data, in_pattern)
!
plan = fftw_plan_dft_1d(length        , in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
call   fftw_execute_dft(plan, in_pattern, out_pattern)
call   fftw_destroy_plan(plan)
out_pattern = out_pattern/sqrt(real(length,kind=PREC_DP)) * yscale ! Scale for H or Q scale
!
CALL mapfftfdtoline(num, dsort, k_data, out_pattern)
!
offxy(iz  ) = offxy(iz-1) +   length
offxy(iz+1) = offxy(iz-1) + 2*length
!
xrange = xmax(idata(kdat)) - xmin(idata(kdat))
xstep  = xrange/REAL(lenc(idata(kdat))-1, kind=PREC_DP)
DO i=1,length
   y(offxy(iz-1)+i) = REAL(k_data(i))*SQRT(xrange*xstep)!    * yscale
   y(offxy(iz  )+i) = IMAG(k_data(i))*SQRT(xrange*xstep)!    * yscale
   x(offxy(iz-1)+i) = x(offxy(idata(1)-1)+i)/(xrange*xstep) * xscale
   x(offxy(iz  )+i) = x(offxy(idata(1)-1)+i)/(xrange*xstep) * xscale
ENDDO
lenc(iz  ) = length
lenc(iz+1) = length
lni (iz)   = .FALSE.
lni (iz+1) = .FALSE.
lh5 (iz)   = .FALSE.
lh5 (iz+1) = .FALSE.
fname(iz  ) = 'Fourier_real'
fname(iz+1) = 'Fourier_imag'
iz = iz + 2
call show_data(iz-2)
call show_data(iz-1)
!
DEALLOCATE(k_data )
DEALLOCATE( in_pattern)
DEALLOCATE(out_pattern)
!
END SUBROUTINE kuplot_do_fft_1D
!
!*****7**************************************************************** 
!
SUBROUTINE kuplot_do_fft_2D(idata, odata, xscale, yscale)
!-
!  Calculate 2D FFT of data set idata
!+
!
USE kuplot_mod
use kuplot_extrema_mod
!
use lib_f90_fftw3
USE map_1dtofield
USE wink_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(2), INTENT(IN)  :: idata     ! Data set number to transform
INTEGER, DIMENSION(2), INTENT(OUT) :: odata     ! Data set number to transform
real(kind=PREC_DP)   , intent(in)  :: xscale
real(kind=PREC_DP)   , intent(in)  :: yscale
!
INTEGER :: i
INTEGER :: kdat
INTEGER :: length             ! Data set length == nx * ny
INTEGER, DIMENSION(3) :: num  ! DATA set dimensions
INTEGER, DIMENSION(3) :: dsort
COMPLEX(KIND=KIND(0.0D0)) , DIMENSION(:),   ALLOCATABLE  :: k_data   ! The Kuplot in data set)
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:), ALLOCATABLE  ::  in_pattern  ! The Curve to be FFT'd
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:), ALLOCATABLE  :: out_pattern  ! The result of   FFT
REAL :: xrange
REAL :: xstep
REAL :: yrange
REAL :: ystep
type(c_ptr) :: plan    ! FFWT3 plan
!
!
!  Determine sequence of array dimensions 
!
IF(idata(1)>0) THEN              ! Real is provided
   kdat   = 1
ELSE                             ! Imag only
   kdat   = 2
ENDIF
num(1) = nx(idata(kdat))
num(2) = ny(idata(kdat))
num(3) = 1
length = num(1)*num(2)
!
dsort(1)      = MAXLOC(num, 1)
num(dsort(1)) = -num(dsort(1))
dsort(2)      = MAXLOC(num, 1)
num(dsort(2)) = -num(dsort(2))
dsort(3)      = MAXLOC(num, 1)
num(dsort(3)) = -num(dsort(3))
num = -num
!
!write(*,*) ' IDATA(1)',     idata(1), idata(2)
!write(*,*) ' LENC(1) ', lenc(idata(1))!, ' :', lenc(lbound(lenc,1):4)
!write(*,*) ' OFFXY   ', offxy(0:2)
!write(*,*) ' OFFZ    ', offz (0:2)
!write(*,*) ' IZ      ', iz
!write(*,*) ' NUM ', num
!write(*,*) ' length ', length
!read(*,*) i
ALLOCATE(k_data (length))
!ALLOCATE(pattern(num(dsort(1)), num(dsort(2)) ))
ALLOCATE( in_pattern(num(dsort(1)), num(dsort(2)) ))
ALLOCATE(out_pattern(num(dsort(1)), num(dsort(2)) ))
IF(idata(1)>0 .and. idata(2)>0) THEN            ! Got real and imag part
   DO i=1,length
      k_data(i) = CMPLX(z(offz(idata(1)-1)+i), z(offz(idata(2)-1)+i))
   ENDDO
ELSEIF(idata(1)>0) THEN                         ! Got real only
   DO i=1,length
      k_data(i) = CMPLX(z(offz(idata(1)-1)+i), 0.0D0)
   ENDDO
ELSEIF(idata(2)>0) THEN                         ! Got imag only
   DO i=1,length
      k_data(i) = CMPLX(0.0D0, z(offz(idata(2)-1)+i))
   ENDDO
ENDIF
!
CALL maptofftfd(num, dsort, k_data, in_pattern)
!
plan = fftw_plan_dft_2d(num(dsort(2)), num(dsort(1)), in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
call   fftw_execute_dft(plan, in_pattern, out_pattern)
call   fftw_destroy_plan(plan)
out_pattern = out_pattern/sqrt(real(num(dsort(1))*num(dsort(2)),kind=PREC_DP)) * yscale ! Scale for H or Q scale
!
CALL mapfftfdtoline(num, dsort, k_data, out_pattern)
!
offxy(iz  ) = offxy(iz-1) +   max(num(1), num(2))
offxy(iz+1) = offxy(iz-1) + 2*max(num(1), num(2))
offz (iz  ) = offz(iz-1) +   length
offz (iz+1) = offz(iz-1) + 2*length
!
xrange = xmax(idata(kdat)) - xmin(idata(kdat))
xstep  = xrange/REAL(nx(idata(kdat))-1,kind=PREC_DP)
yrange = ymax(idata(kdat)) - ymin(idata(kdat))
ystep  = yrange/REAL(ny(idata(kdat))-1, kind=PREC_DP)
!
DO i=1,length
   z(offz(iz-1)+i) = REAL(k_data(i))*SQRT(xrange*xstep*yrange*ystep)
   z(offz(iz  )+i) = IMAG(k_data(i))*SQRT(xrange*xstep*yrange*ystep)
ENDDO
DO i =1, nx(idata(kdat))
   x(offxy(iz-1) +i) = x(offxy(idata(kdat)-1)+i)/(xrange*xstep) * xscale
   x(offxy(iz  ) +i) = x(offxy(idata(kdat)-1)+i)/(xrange*xstep) * xscale
ENDDO
DO i =1, ny(idata(kdat))
   y(offxy(iz-1) +i) = y(offxy(idata(kdat)-1)+i)/(yrange*ystep) * xscale
   y(offxy(iz  ) +i) = y(offxy(idata(kdat)-1)+i)/(yrange*ystep) * xscale
ENDDO
fform(iz  ) = 'NI'
fform(iz+1) = 'NI'
nx  (iz)   = nx(idata(kdat))
ny  (iz)   = ny(idata(kdat))
nx  (iz+1) = nx(idata(kdat))
ny  (iz+1) = ny(idata(kdat))
lni (iz)   = .TRUE.
lni (iz+1) = .TRUE.
lh5 (iz)   = .FALSE.
lh5 (iz+1) = .FALSE.
lenc(iz  ) = MAX(nx  (iz), ny(iz))
lenc(iz+1) = MAX(nx  (iz), ny(iz))
fname(iz  ) = 'Fourier_real'
fname(iz+1) = 'Fourier_imag'
CALL get_extrema_xy (x, iz, nx (iz), xmin, xmax)
CALL get_extrema_xy (y, iz, ny (iz), ymin, ymax)
CALL get_extrema_z  (z, iz, nx (iz), ny (iz), zmin, zmax)
CALL get_extrema_xy (x, iz+1, nx (iz+1), xmin, xmax)
CALL get_extrema_xy (y, iz+1, ny (iz+1), ymin, ymax)
CALL get_extrema_z  (z, iz+1, nx (iz+1), ny (iz+1), zmin, zmax)
iz = iz + 2
!write(*,*) ' OFFXY   ', offxy(0:4)
!write(*,*) ' OFFZ    ', offz (0:4)
!write(*,*) ' IZ      ', iz
!
DEALLOCATE(k_data)
DEALLOCATE( in_pattern)
DEALLOCATE(out_pattern)
!
END SUBROUTINE kuplot_do_fft_2D
!
!*****7**************************************************************** 
!
!
!*****7**************************************************************** 
!
SUBROUTINE do_glat (zeile, lp, lsmooth) 
!+                                                                      
!     Smooth data set                                                   
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
USE precision_mod
USE string_convert_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw, maxsm 
PARAMETER (maxw = 4) 
PARAMETER (maxsm = 501) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lp
LOGICAL         , intent(in)    :: lsmooth 
!                                                                       
CHARACTER(len=2)                             :: dir 
CHARACTER(LEN=PREC_STRING), dimension(MAXW)  :: cpara (maxw) 
INTEGER                   , dimension(MAXW)  :: lpara (maxw) 
REAL(KIND=PREC_DP)        , dimension(MAXW)  :: werte (maxw) 
REAL(kind=PREC_SP)        , dimension(MAXSM) :: cc (maxsm) 
INTEGER :: ianz, ik, ip, im, il 
LOGICAL :: gx, gy 
!                                                                       
IF (lsmooth) then 
   DO ip = 1, maxsm 
      cc (ip) = 0.0 
   ENDDO 
ELSE 
   DO ip = 1, maxsm 
      cc (ip) = 1.0 
   ENDDO 
ENDIF 
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.eq.2.or.ianz.eq.3.or.ianz.eq.4) then 
   IF (lsmooth.and.ianz.eq.4) then 
      CALL do_cap (cpara (4) ) 
      gx = (cpara (4) (1:1) .eq.'X') 
      gy = (cpara (4) (1:1) .eq.'Y') 
      cpara (4) = '0.0' 
   ELSEIF (.not.lsmooth.and.ianz.eq.3) then 
      CALL do_cap (cpara (3) ) 
      gx = (cpara (3) (1:1) .eq.'X') 
      gy = (cpara (3) (1:1) .eq.'Y') 
      cpara (3) = '0.0' 
   ELSE 
      gx = .true. 
      gy = .true. 
   ENDIF 
!                                                                       
   CALL ber_params (ianz, cpara, lpara, werte, maxw) 
   ik = nint (werte (1) ) 
   ip = nint (werte (2) ) 
   IF (ik.lt.iz.and.ik.gt.0) then 
!                                                                       
      IF (gx) dir = ' x' 
      IF (gy) dir = ' y' 
      IF (gx.and.gy) dir = 'xy' 
!                                                                       
      IF (lsmooth) then 
         il = 0 
         im = 2 
         IF (ianz.gt.2) im = nint (werte (3) ) 
         CALL SAVGOL (cc, ip, ip / 2, ip / 2, il, im) 
         IF (ier_num.ne.0) return 
      ENDIF 
!                                                                       
      IF (lni (ik) ) then 
         IF (lsmooth) then 
            WRITE (output_io, 2000) ik, ip, im 
         ELSE 
            WRITE (output_io, 2010) ik, ip 
         ENDIF 
         WRITE (output_io, 2100) dir 
         CALL do_glatt_z(ik, ip, nx (ik), ny (ik), gx, gy, cc, maxsm, lsmooth)
      ELSE 
         IF (lsmooth) then 
            WRITE (output_io, 2000) ik, ip, im 
         ELSE 
            WRITE (output_io, 2010) ik, ip 
         ENDIF 
         CALL do_glatt_y (ik, ip, lenc(ik), cc, maxsm, lsmooth) 
      ENDIF 
   ELSE 
      ier_num = - 4 
      ier_typ = ER_APPL 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
 2000 FORMAT (' Smoothing data set ',i3,' via Savitzky-Golay :',/,      &
     &        '   Number of points           : ',i3,/,                  &
     &        '   Order of smoothing polynom : ',i3)                    
 2010 FORMAT (' Smoothing data set ',i3,' via sliding average :',/,     &
     &        '   Number of points           : ',i3)                    
 2100 FORMAT ('   Direction                  : ',1x,a2) 
END SUBROUTINE do_glat                        
!
!*****7*****************************************************************
!
SUBROUTINE do_inte (zeile, lp) 
!+                                                                      
!     Compute integral for 2d/3d data sets                              
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE param_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_para_mod
!
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 5) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lp
!
CHARACTER(LEN=PREC_STRING), dimension(MAXW) :: cpara (maxw) 
INTEGER                   , dimension(MAXW) :: lpara (maxw)
INTEGER ::ianz, ipkt, ip, ik, ix, iy 
REAL(kind=PREC_DP) :: rint, drint, rdx, rdy, xa, ya, xs 
REAL(kind=PREC_DP) :: wx1, wx2, wy1, wy2 
REAL(KIND=PREC_DP)        , dimension(MAXW) :: werte (maxw) 
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ianz.lt.1) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
IF (ier_num.ne.0) return 
CALL skalieren 
ik = nint (werte (1) ) 
wx1 = ex (iwin, iframe, 1) 
wx2 = ex (iwin, iframe, 2) 
wy1 = ey (iwin, iframe, 1) 
wy2 = ey (iwin, iframe, 2) 
!                                                                       
IF (ik.gt.0.and.ik.lt.iz) then 
   IF (ianz.eq.1) then 
      CONTINUE 
   ELSEIF (ianz.eq.3) then 
      wx1 = werte (2) 
      wx2 = werte (3) 
   ELSEIF (ianz.eq.5.and.lni (ik) ) then 
      wx1 = werte (2) 
      wx2 = werte (3) 
      wy1 = werte (4) 
      wy2 = werte (5) 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
ELSE 
   ier_num = - 4 
   ier_typ = ER_APPL 
ENDIF 
!                                                                       
IF (ier_num.ne.0) return 
!                                                                       
rint = 0.0 
drint = 0.0 
ipkt = 0 
!                                                                       
!------ 2D integration                                                  
!                                                                       
IF (lni (ik) ) then 
   DO ix = 1, nx (ik) 
      DO iy = 1, ny (ik) 
         IF (x (offxy (ik - 1) + ix) .ge.wx1.and.y (offxy (ik - 1)      &
         + iy) .ge.wy1.and.x (offxy (ik - 1) + ix) .le.wx2.and.y (offxy &
         (ik - 1) + iy) .le.wy2) then                                   
            zwert = z (offz (ik - 1) + (ix - 1) * ny (ik) + iy) 
            IF (zwert.ne. - 9999.0) then 
               rint = rint + zwert 
               ipkt = ipkt + 1 
            ENDIF 
         ENDIF 
      ENDDO 
   ENDDO 
   rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) ) 
   rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) ) 
   WRITE (output_io, 2000) ik, wx1, wx2, wy1, wy2, rint * rdx * rdy, drint, ipkt
   IF (maxpar_res.lt.7) then 
      ier_num = - 24 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
   res_para (0) = 7 
   res_para (1) = rint * rdx * rdy 
   res_para (2) = drint 
   res_para (3) = REAL(ipkt) 
   res_para (4) = wx1 
   res_para (5) = wx2 
   res_para (6) = wy1 
   res_para (7) = wy2 
!                                                                       
!------ 1D integration                                                  
!                                                                       
ELSE 
   DO ix = 2, lenc(ik) 
      ip = offxy (ik - 1) + ix 
      xa = (x (ip) + x (ip - 1) ) / 2.0 
      ya = (y (ip) + y (ip - 1) ) / 2.0 
      xs = abs (x (ip) - x (ip - 1) ) 
!                                                                       
      IF (xa.ge.wx1.and.xa.le.wx2) then 
         rint = rint + xs * ya 
         drint = drint + xs**2 * (dy(ip) **2 + dy(ip - 1) **2) / 4.0 
         ipkt = ipkt + 1 
      ENDIF 
   ENDDO 
   WRITE (output_io, 2100) ik, wx1, wx2, rint, sqrt (drint), ipkt
   IF (maxpar_res.lt.5) then 
      ier_num = - 24 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
   res_para (0) = 5 
   res_para (1) = rint 
   res_para (2) = sqrt (drint) 
   res_para (3) = REAL(ipkt) 
   res_para (4) = wx1 
   res_para (5) = wx2 
ENDIF 
WRITE (output_io, * ) 
!                                                                       
 2000 FORMAT    (' Integration result for data set ',i3,' :',/          &
     &                  3x,'x-range    : ',g12.4,' to ',g12.4,/         &
     &                  3x,'y-range    : ',g12.4,' to ',g12.4,/         &
     &         3x,'Integral   : ',g12.4,' +- ',g12.4,'  (',i7,' pkt)')  
 2100 FORMAT    (' Integration result for data set ',i3,' :',/          &
     &                  3x,'x-range    : ',g12.4,' to ',g12.4,/         &
     &         3x,'Integral   : ',g12.4,' +- ',g12.4,'  (',i7,' pkt)')  
!                                                                       
END SUBROUTINE do_inte                        
!
!**7******************************************************************* 
!
SUBROUTINE do_mean (zeile, lp, lout) 
!+                                                                      
!     calculates mean and standard deviation                            
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE param_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_para_mod
!
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 5) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lp
LOGICAL         , intent(in)    :: lout 
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
INTEGER :: lpara (maxw)
INTEGER :: ianz, ik, ix, iy, i 
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL(kind=PREC_DP) :: wx1, wx2, wy1, wy2 
REAL(kind=PREC_DP) :: sum (3), sum2 (3), wert (3), ave (3), sig (3), n 
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ianz.lt.1) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
IF (ier_num.ne.0) return 
!                                                                       
ik = nint (werte (1) ) 
IF (ik.ge.iz.or.lenc(ik) .eq.0) then 
   ier_num = - 4 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
CALL skalieren 
IF (ianz.eq.1) then 
   wx1 = ex (iwin, iframe, 1) 
   wx2 = ex (iwin, iframe, 2) 
   wy1 = ey (iwin, iframe, 1) 
   wy2 = ey (iwin, iframe, 2) 
ELSEIF (ianz.eq.3) then 
   wx1 = werte (2) 
   wx2 = werte (3) 
   wy1 = ey (iwin, iframe, 1) 
   wy2 = ey (iwin, iframe, 2) 
ELSEIF (ianz.eq.5.and.lni (ik) ) then 
   wx1 = werte (2) 
   wx2 = werte (3) 
   wy1 = werte (4) 
   wy2 = werte (5) 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
!DO i = 1, 3 
!   sum (i) = 0.0 
!   sum2 (i) = 0.0 
!ENDDO 
sum  = 0.0D0
sum2 = 0.0D0
n = 0 
!                                                                       
IF (lni (ik) ) then 
   DO ix = 1, nx (ik) 
      DO iy = 1, ny (ik) 
         IF (x (offxy (ik - 1) + ix) .ge.wx1.and.y (offxy (ik - 1)      &
         + iy) .ge.wy1.and.x (offxy (ik - 1) + ix) .le.wx2.and.y (offxy &
         (ik - 1) + iy) .le.wy2) then                                   
            wert (1) = x (offxy (ik - 1) + ix) 
            wert (2) = y (offxy (ik - 1) + iy) 
            wert (3) = z (offz (ik - 1) + (ix - 1) * ny (ik) + iy) 
            IF (werte (3) .ne. - 9999.0) then 
               DO i = 1, 3 
               sum (i) = sum (i) + wert (i) 
               sum2 (i) = sum2 (i) + wert (i) **2 
               ENDDO 
               n = n + 1 
            ENDIF 
         ENDIF 
      ENDDO 
   ENDDO 
   IF (n.gt.0) then 
      DO i = 1, 3 
         ave (i) = sum (i) / n 
         sig (i) = sqrt (sum2 (i) / n - (sum (i) / n) **2) 
      ENDDO 
      IF (lout) write (output_io, 2000) ik, wx1, wx2, wy1, wy2,   &
            int (n), (ave (i), sig (i), i = 1, 3)                       
      IF (maxpar_res.lt.11) then 
         ier_num = - 24 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      res_para (0) = 11 
      res_para (1) = ave (1) 
      res_para (2) = ave (2) 
      res_para (3) = ave (3) 
      res_para (4) = sig (1) 
      res_para (5) = sig (2) 
      res_para (6) = sig (3) 
      res_para (7) = int (n) 
      res_para (8) = wx1 
      res_para (9) = wx2 
      res_para (10) = wy1 
      res_para (11) = wy2 
   ELSE 
      ier_num = - 20 
      ier_typ = ER_APPL 
   ENDIF 
ELSE 
   DO ix = 1, lenc(ik) 
      IF(x(offxy(ik - 1) + ix) .ge. wx1 .and. x(offxy(ik - 1) + ix) .le.wx2) then
         wert (1) = x (offxy (ik - 1) + ix) 
         wert (2) = y (offxy (ik - 1) + ix) 
         DO i = 1, 2 
            sum (i) = sum (i) + wert (i) 
            sum2 (i) = sum2 (i) + wert (i) **2 
         ENDDO 
         n = n + 1 
      ENDIF 
   ENDDO 
   IF (n.gt.0) then 
      DO i = 1, 2 
         ave (i) = sum (i) / n 
         sig (i) = sqrt (sum2 (i) / n - (sum (i) / n) **2) 
      ENDDO 
      IF (lout) write (output_io, 2100) ik, wx1, wx2, int (n),    &
            (ave (i), sig (i), i = 1, 2)                                
      IF (maxpar_res.lt.7) then 
         ier_num = - 24 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      res_para (0) = 7 
      res_para (1) = ave (1) 
      res_para (2) = ave (2) 
      res_para (3) = sig (1) 
      res_para (4) = sig (2) 
      res_para (5) = int (n) 
      res_para (6) = wx1 
      res_para (7) = wx2 
   ELSE 
      ier_num = - 20 
      ier_typ = ER_APPL 
   ENDIF 
ENDIF 
IF (lout) write (output_io, * ) 
!                                                                       
 2000 FORMAT    (' Average and standard deviation data set ',i3,' :',/  &
     &                  3x,'x-range   : ',g12.4,' to ',g12.4,/          &
     &                  3x,'y-range   : ',g12.4,' to ',g12.4,/          &
     &                  3x,'# points  : ',i7//                          &
     &         3x,'Average x : ',g12.4,' sigma :  ',g12.4,/             &
     &         3x,'Average y : ',g12.4,' sigma :  ',g12.4,/             &
     &         3x,'Average z : ',g12.4,' sigma :  ',g12.4)              
 2100 FORMAT    (' Average and standard deviation data set ',i3,' :',/  &
     &                  3x,'x-range   : ',g12.4,' to ',g12.4,/          &
     &                  3x,'# points  : ',i7//                          &
     &         3x,'Average x : ',g12.4,' sigma :  ',g12.4,/             &
     &         3x,'Average y : ',g12.4,' sigma :  ',g12.4)              
!                                                                       
END SUBROUTINE do_mean                        
!
!**7******************************************************************* 
!
SUBROUTINE do_smax (zeile, lp) 
!+                                                                      
!     Search for maxima in data set                                     
!-                                                                      
USE ber_params_mod
USE build_name_mod
USE errlist_mod 
USE get_params_mod
USE param_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
USE lib_length
USE precision_mod
USE support_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxmax, maxw 
PARAMETER (maxmax = 50) 
PARAMETER (maxw = 10) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
CHARACTER(len=80)          :: cdummy 
CHARACTER(len=6)           :: cout 
REAL(KIND=PREC_DP), dimension(MAXW)   :: werte (maxw) 
REAL(kind=PREC_DP), dimension(MAXMAX) :: wmax (maxmax) 
INTEGER :: lpara (maxw)
INTEGER :: ixm (maxmax), iym (maxmax) 
INTEGER :: ianz, ifil, im, ima, ik 
!                                                                       
!                                                                       
cout = 'maxima' 
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.ge.2) then 
   CALL ber_params (2, cpara, lpara, werte, maxw) 
   IF (ier_num.ne.0) return 
   ik = nint (werte (1) ) 
   ifen = nint (werte (2) ) 
   IF (ik.gt.0.and.ik.lt.iz) then 
      IF (ianz.gt.2) then 
         CALL del_params (2, ianz, cpara, lpara, maxw) 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) return 
         ifil = 66 
         cdummy = cpara (1) (1:MIN(80,LEN_TRIM(cpara(1))))
         CALL oeffne (ifil, cdummy, 'unknown') 
         IF (ier_num.ne.0) return 
         WRITE (output_io, 900) cdummy (1:len_str (cdummy) ) 
      ELSE 
         ifil = output_io 
      ENDIF 
!                                                                       
!------ --- 3D data set                                                 
!                                                                       
      IF (lni (ik) ) then 
         CALL do_fmax_z (ik, wmax, ixm, iym, maxmax, ima) 
         IF (ier_num.eq.0) then 
            WRITE (ifil, 1000) cout, ik, ifen 
            IF (maxpar_res.lt.2 * ima) then 
               ier_num = - 24 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
            res_para (0) = 2.0 * ima 
            DO im = 1, ima 
               WRITE (ifil, 1100) im, x (offxy (ik - 1) + ixm (im) ),&
               y (offxy (ik - 1) + iym (im) ), wmax (im)             
               res_para (2 * im - 1) = REAL(ixm (im) ) 
               res_para (2 * im) = REAL(iym (im) ) 
            ENDDO 
         ENDIF 
!                                                                       
!------ --- 2D data set                                                 
!                                                                       
      ELSE 
         CALL do_fmax_xy (ik, wmax, ixm, maxmax, ima) 
         IF (ier_num.eq.0) then 
            WRITE (ifil, 2000) cout, ik, ifen 
            IF (maxpar_res.lt.ima) then 
               ier_num = - 24 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
            res_para (0) = ima 
            DO im = 1, ima 
               WRITE (ifil, 2100) im, x (offxy (ik - 1) + ixm (im) ),&
               wmax (im)                                             
               res_para (im) = REAL(ixm (im) ) 
            ENDDO 
         ENDIF 
      ENDIF 
      IF (ifil.ne.6) close (ifil) 
   ELSE 
      ier_num = - 4 
      ier_typ = ER_APPL 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
  900 FORMAT (' ------ > Saving results to file ',a,' ..') 
 1000 FORMAT (' Found ',a,' data set ',i3,' (ifen = ',i3,') :',//,      &
     &        3x,'No.        pos. x       pos. y         value'/        &
     &        3x,'--------------------------------------------------')  
 1100 FORMAT (3x,i3,2x,f11.3,2x,f11.3,2x,f18.3) 
 2000 FORMAT (' Found ',a,' data set ',i3,' (ifen = ',i3,') :',//,      &
     &        3x,'No.        pos. x         value'/                     &
     &        3x,'-------------------------------------')               
 2100 FORMAT (3x,i3,2x,f11.3,2x,f18.3) 
!
END SUBROUTINE do_smax                        
!
!**7*****************************************************************   
!     SUPPORT ROUTINES                                                  
!**7*****************************************************************   
!
SUBROUTINE do_glatt_z (ik, np, nxx, nyy, gx, gy, cc, mc, lsm) 
!+                                                                      
!     glaetten des nipl-2d-arrays                                       
!-                                                                      
USE kuplot_config 
USE kuplot_mod 
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, intent(in) :: ik
INTEGER, intent(in) :: np
INTEGER, intent(in) :: nxx
INTEGER, intent(in) :: nyy
LOGICAL, intent(in) :: gx
LOGICAL, intent(in) :: gy
INTEGER, intent(in) :: MC
REAL(kind=PREC_SP), intent(in) :: cc (mc) 
LOGICAL, intent(in) :: lsm 
!                                                                       
!REAL(kind=PREC_DP), DIMENSION(MAXARRAY) :: a (maxarray)
REAL(kind=PREC_DP), DIMENSION(:), allocatable :: a !(maxarray)
REAL(kind=PREC_DP) :: xnw, s 
!     INTEGER ik, np, nxx, nyy 
INTEGER :: nd, na, ne, ip, jp, im, jm, i, j, k 
!                                                                       
nd = (np - 1) / 2 
na = nd+1 
!                                                                       
!------ glaettung in x-richtung                                         
!                                                                       
!     IF (.not.gx) goto 222 
if_gx: if(gx) then
   ne = nxx - nd 
   allocate(a(na:ne))
   DO j = 1, nyy 
      DO i = na, ne 
         xnw = np 
         im = i - nd 
         ip = i + nd 
         s = 0.0 
         DO k = im, ip 
            IF (z (offz (ik - 1) + (k - 1) * ny (ik) + j) .ne. - 9999.0.and.z &
               (offz (ik - 1) + (k - 1) * ny (ik) + j) .ne.0.0) then             
               s = s + cc (abs (i - k) + 1) * z (offz (ik - 1) + (k - 1)      &
               * ny (ik) + j)                                                 
            ELSE 
               xnw = xnw - 1 
            ENDIF 
         ENDDO 
         IF (xnw.gt.0) then 
            IF (lsm) then 
               a (i) = s 
            ELSE 
               a (i) = s / xnw 
            ENDIF 
         ELSE 
            a (i) = - 9999.0 
         ENDIF 
      ENDDO 
      DO i = na, ne 
         z (offz (ik - 1) + (i - 1) * ny (ik) + j) = a (i) 
      ENDDO 
   ENDDO 
   deallocate(a)
endif if_gx
!                                                                       
!------ glaettung in y-richtung                                         
!                                                                       
!  222 CONTINUE 
!     IF (.not.gy) goto 333 
if_gy: if(gy) then
   ne = nyy - nd 
   allocate(a(na:ne))
   DO i = 1, nxx 
      DO j = na, ne 
         xnw = np 
         jm = j - nd 
         jp = j + nd 
         s = 0.0 
         DO k = jm, jp 
            IF(z(offz(ik - 1) + (i - 1) * ny(ik) + k) .ne. - 9999.0 .and.     & 
               z(offz(ik - 1) + (i - 1) * ny (ik) + k) .ne.0.0            ) then             
               s = s + cc(abs(j - k) + 1) * z(offz(ik - 1) + (i - 1) * ny (ik) + k)
            ELSE 
               xnw = xnw - 1 
            ENDIF 
         ENDDO 
         IF (xnw.gt.0) then 
            IF (lsm) then 
               a (j) = s 
            ELSE 
               a (j) = s / xnw 
            ENDIF 
         ELSE 
            a (j) = - 9999.0 
         ENDIF 
      ENDDO 
      DO j = na, ne 
         z (offz (ik - 1) + (i - 1) * ny (ik) + j) = a (j) 
      ENDDO 
   ENDDO 
   deallocate(a)
endif if_gy
!                                                                       
! 333 CONTINUE 
END SUBROUTINE do_glatt_z                     
!
!**7*****************************************************************   
!
SUBROUTINE do_glatt_y (ik, ip, llen, cc, mc, lsm) 
!+                                                                      
!     glaetten der kurve ik ueber ip punkte                             
!-                                                                      
USE kuplot_config 
USE kuplot_mod 
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, intent(in) :: ik
INTEGER, intent(in) :: ip
INTEGER, intent(in) :: llen
INTEGER, intent(in) :: MC
REAL(kind=PREC_SP), dimension(MC), intent(in)  :: cc (mc) 
LOGICAL, intent(in) :: lsm 
!
!REAL(kind=PREC_DP) :: ay (maxarray), ady (maxarray)! , cc (mc) 
REAL(kind=PREC_DP), dimension(:), allocatable :: ay  !(nd+1:ne)
REAL(kind=PREC_DP), dimension(:), allocatable :: ady !(nd+1:ne)! , cc (mc) 
REAL(kind=PREC_DP) :: s, ds, xnw 
INTEGER :: na, ne, nd, ic, i, k, kk 
!                                                                       
nd = (ip - 1) / 2 
na = nd+1 
ne = llen - nd 
allocate(ay(na:ne))
allocate(ady(na:ne))
!                                                                       
DO i = na, ne 
   xnw = ip 
   s = 0.0 
   ds = 0.0 
!                                                                       
   DO k = - nd, 0 
      kk = k + i 
      ic = 1 - k 
      s = s + cc (ic) * y (offxy (ik - 1) + kk) 
      ds = ds + cc (ic) **2 * dy (offxy (ik - 1) + kk) **2 
   ENDDO 
!                                                                       
   DO k = 1, nd 
      kk = k + i 
      ic = 2 * nd+2 - k 
      s = s + cc (ic) * y (offxy (ik - 1) + kk) 
      ds = ds + cc (ic) **2 * dy (offxy (ik - 1) + kk) **2 
   ENDDO 
!                                                                       
   IF (lsm) then 
      ay (i) = s 
      ady (i) = sqrt (ds) 
   ELSE 
      ay (i) = s / REAL(ip) 
      ady (i) = sqrt (ds / REAL(ip) ) 
   ENDIF 
ENDDO 
!                                                                       
DO i = na, ne 
   y (offxy (ik - 1) + i) = ay (i) 
   dy (offxy (ik - 1) + i) = ady (i) 
ENDDO 
!
deallocate(ay) 
deallocate(ady)
!                                                                       
END SUBROUTINE do_glatt_y                     
!
!**7*****************************************************************   
!
SUBROUTINE extract_subarray(xf, yf, zf, xmit, ymit, nsize, ik, ie)
!+                                                                      
!     schneiden eines subarrays fuer interpolationsroutinen             
!-                                                                      
USE kuplot_config 
USE kuplot_mod 
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER                                    , intent(IN) :: nsize 
real(kind=PREC_SP), dimension(NSIZE)       , intent(OUT) :: xf
real(kind=PREC_SP), dimension(NSIZE)       , intent(OUT) :: yf
real(kind=PREC_SP), dimension(NSIZE, NSIZE), intent(OUT) :: zf
!     REAL xf (nsize), yf (nsize), zf (nsize, nsize) 
REAL(kind=PREC_SP)                         , intent(in) :: xmit
REAL(kind=PREC_SP)                         , intent(in) :: ymit 
integer                                    , intent(in)  :: ik
integer                                    , intent(out) :: ie
!
INTEGER :: ip, ii, jj, ixx, iyy 
!                                                                       
ie = 0 
ip = nsize / 2 + 1 
ixx = nint((xmit - xmin(ik)) / (xmax(ik) - xmin(ik)) * REAL(nx(ik) - 1)) + 1
iyy = nint((ymit - ymin(ik)) / (ymax(ik) - ymin(ik)) * REAL(ny(ik) - 1)) + 1
IF(ixx.lt.ip.or. (ixx + ip) .gt.nx (ik) .or.iyy.lt.ip.or. (iyy + ip) .gt.ny (ik) ) then
   ie = - 1 
   RETURN 
ENDIF 
DO ii = 1, nsize 
   DO jj = 1, nsize 
      xf (ii) = x (offxy (ik - 1) + ii - ip + ixx) 
      yf (jj) = y (offxy (ik - 1) + jj - ip + iyy) 
      zf (ii, jj) = z (offz (ik - 1) + (ii - ip + ixx) * ny (ik)        &
      + jj - ip + iyy)                                                  
      IF (zf (ii, jj) .eq. - 9999.) then 
         ie = - 1 
         RETURN 
      ENDIF 
   ENDDO 
ENDDO 
!                                                                       
END SUBROUTINE extract_subarray               
!
!********************************************************************   
!
SUBROUTINE polin2 (x1a, x2a, ya, m, n, x1, x2, y, dy,ier) 
!
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
!
integer                , intent(in) :: m
integer                , intent(in) :: n
REAL   , dimension(m)  , intent(in) :: x1a
REAL   , dimension(n)  , intent(in) :: x2a
REAL   , dimension(m,n), intent(in) :: ya
REAL                   , intent(in) :: x1
REAL                   , intent(in) :: x2
REAL                   , intent(out) :: y
REAL                   , intent(out) :: dy
integer                , intent(out) :: ier
!
integer, PARAMETER :: nmax = 50
integer, PARAMETER :: mmax = 50
!      DIMENSION x1a (m), x2a (n), ya (m, n), 
real :: yntmp (nmax), ymtmp (mmax) 
!      INTEGER :: ier
      ier = 0
      DO 12 j = 1, m 
         DO 11 k = 1, n 
            yntmp (k) = ya (j, k) 
   11    END DO 
         CALL polint (x2a, yntmp, n, x2, ymtmp (j), dy, ier) 
         IF(ier/=0) RETURN
   12 END DO 
      CALL polint (x1a, ymtmp, m, x1, y, dy, ier) 
      RETURN 
!
END SUBROUTINE polin2                         
!
!********************************************************************   
!
SUBROUTINE polint (xa, ya, n, x, y, dy, ier) 
!
use precision_mod
!
!IMPLICIT integer(i-n)
!IMPLICIT REAL (a-h, o-z)
implicit none
!
integer                          , intent(in) :: N
real(kind=PREC_SP) , dimension(N), intent(in) :: xa
real(kind=PREC_SP) , dimension(N), intent(in) :: ya
real(kind=PREC_SP)               , intent(in)  :: x
real(kind=PREC_SP)               , intent(out) :: y
real(kind=PREC_SP)               , intent(out) :: dy
INTEGER                          , intent(out) :: ier

integer, PARAMETER :: nmax = 50 
!
integer ::  i
integer ::  m
integer ::  ns
real(kind=PREC_DP), dimension(NMAX) :: c
real(kind=PREC_DP), dimension(NMAX) :: d
real(kind=PREC_DP) :: dif
real(kind=PREC_DP) :: dift
real(kind=PREC_DP) :: den
real(kind=PREC_DP) :: ho
real(kind=PREC_DP) :: hp
real(kind=PREC_DP) :: w
!     DIMENSION xa (n), ya (n), c (nmax), d (nmax) 
!
ier = 0
ns = 1 
dif = abs (x - xa (1) ) 
!
DO i = 1, n 
   dift = abs (x - xa (i) ) 
   IF (dift.lt.dif) then 
      ns = i 
      dif = dift 
   ENDIF 
   c (i) = ya (i) 
   d (i) = ya (i) 
END DO 
y = ya (ns) 
ns = ns - 1 
!
DO m = 1, n - 1 
   DO i = 1, n - m 
      ho = xa (i) - x 
      hp = xa (i + m) - x 
      w = c (i + 1) - d (i) 
      den = ho - hp 
      IF (den.eq.0.) THEN
         ier = -60
         RETURN
      ENDIF
      den = w / den 
      d (i) = hp * den 
      c (i) = ho * den 
   END DO 
!
   IF (2 * ns.lt.n - m) then 
      dy = c (ns + 1) 
   ELSE 
      dy = d (ns) 
      ns = ns - 1 
   ENDIF 
   y = y + dy 
END DO 
!     RETURN 
!
END SUBROUTINE polint                         
!                                                                       
!*******************************************************************************
!
SUBROUTINE spline_old (x, y, n, yp1, ypn, y2) 
!
USE kuplot_config 
use precision_mod
!
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
!     PARAMETER (nmax = maxarray) 
      DIMENSION x (n), y (n), y2 (n) 
!dimension u (nmax) 
dimension u (n   ) 
      IF (yp1.gt..99e30) then 
         y2 (1) = 0. 
         u (1) = 0. 
      ELSE 
         y2 (1) = - 0.5 
         u (1) = (3. / (x (2) - x (1) ) ) * ( (y (2) - y (1) ) /        &
         (x (2) - x (1) ) - yp1)                                        
      ENDIF 
      DO 11 i = 2, n - 1 
         sig = (x (i) - x (i - 1) ) / (x (i + 1) - x (i - 1) ) 
         p = sig * y2 (i - 1) + 2. 
         y2 (i) = (sig - 1.) / p 
         u (i) = (6. * ( (y (i + 1) - y (i) ) / (x (i + 1) - x (i) )    &
         - (y (i) - y (i - 1) ) / (x (i) - x (i - 1) ) ) / (x (i + 1)   &
         - x (i - 1) ) - sig * u (i - 1) ) / p                          
   11 END DO 
      IF (ypn.gt..99e30) then 
         qn = 0. 
         un = 0. 
      ELSE 
         qn = 0.5 
         un = (3. / (x (n) - x (n - 1) ) ) * (ypn - (y (n) - y (n - 1) )&
         / (x (n) - x (n - 1) ) )                                       
      ENDIF 
      y2 (n) = (un - qn * u (n - 1) ) / (qn * y2 (n - 1) + 1.) 
      DO 12 k = n - 1, 1, - 1 
         y2 (k) = y2 (k) * y2 (k + 1) + u (k) 
   12 END DO 
      RETURN 
      END SUBROUTINE spline_old                         
!                                                                       
      SUBROUTINE splint_old (xa, ya, y2a, n, x, y, ier) 
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
      DIMENSION xa (n), ya (n), y2a (n) 
      INTEGER :: ier
      klo = 1 
      khi = n 
    1 IF (khi - klo.gt.1) then 
         k = (khi + klo) / 2 
         IF (xa (k) .gt.x) then 
            khi = k 
         ELSE 
            klo = k 
         ENDIF 
         GOTO 1 
      ENDIF 
      h = xa (khi) - xa (klo) 
      IF (h.eq.0.) THEN
         ier = -60
         RETURN
      ENDIF
      a = (xa (khi) - x) / h 
      b = (x - xa (klo) ) / h 
      y = a * ya (klo) + b * ya (khi) + ( (a**3 - a) * y2a (klo)        &
      + (b**3 - b) * y2a (khi) ) * (h**2) / 6.                          
      RETURN 
      END SUBROUTINE splint_old                         
!
!*******************************************************************************
!
SUBROUTINE SAVGOL (c, np, nl, nr, ld, m) 
!
USE errlist_mod 
use precision_mod
!
!IMPLICIT integer(i-n)
!IMPLICIT REAL (a-h, o-z)
implicit none
!
!
integer, intent(in) :: np
real(kind=PREC_SP), dimension(np), intent(out) ::c
integer, intent(in) :: nl
integer, intent(in) :: nr
integer, intent(in) :: ld
integer, intent(in) :: m 
!      INTEGER ld, m, nl, np, nr, MMAX 
!     REAL c (np) 
integer, PARAMETER :: MMAX = 6 
!
INTEGER :: imj, ipj, j, k, kk, mm, indx (MMAX + 1) 
REAL(kind=PREC_DP) :: d, fac, sum, a (MMAX + 1, MMAX + 1), b (MMAX + 1) 
!                                                                       
IF (np.lt.nl + nr + 1 .or. nl.lt.0 .or. nr.lt.0 .or. ld.gt.m .or. & 
    m.gt.MMAX .or. nl + nr.lt.m  .or.  np.le.3) then                                          
   ier_num = - 39 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
DO ipj = 0, 2 * m 
   sum = 0. 
   IF (ipj.eq.0) sum = 1. 
   DO k = 1, nr 
      sum = sum + REAL(k, kind=PREC_DP) **ipj 
   ENDDO 
   DO k = 1, nl 
      sum = sum + REAL(-k, kind=PREC_DP) **ipj 
   ENDDO 
   mm = min (ipj, 2 * m - ipj) 
   DO imj = - mm, mm, 2 
      a (1 + (ipj + imj) / 2, 1 + (ipj - imj) / 2) = sum 
   ENDDO 
ENDDO 
CALL ludcmp (a, m + 1, MMAX + 1, indx, d, ier_num)
IF(ier_num /= 0) THEN
   ier_typ = ER_NONE
   RETURN
ENDIF 
DO j = 1, m + 1 
   b (j) = 0. 
ENDDO 
b (ld+1) = 1. 
CALL lubksb (a, m + 1, MMAX + 1, indx, b) 
DO kk = 1, np 
   c (kk) = 0. 
ENDDO 
DO k = - nl, nr 
   sum = b (1) 
   fac = 1. 
   DO mm = 1, m 
      fac = fac * k 
      sum = sum + b (mm + 1) * fac 
   ENDDO 
   kk = mod (np - k, np) + 1 
   c (kk) = sum 
ENDDO 
!      RETURN 
!
END SUBROUTINE SAVGOL                         
!
!*******************************************************************************
!
SUBROUTINE LUBKSB (A, N, NP, INDX, B) 
!
use precision_mod
!
implicit none
!IMPLICIT integer(i-n)
!IMPLICIT REAL (a-h, o-z)
!      DIMENSION A (NP, NP), INDX (N), B (N) 
!
integer                              , intent(in)    :: N
integer                              , intent(in)    :: NP
real(kind=PREC_DP), dimension(NP, NP), intent(inout) :: a
integer           , dimension(N)     , intent(out)   :: indx
real(kind=PREC_DP), dimension(NP)    , intent(inout) :: b
!
integer :: i
integer :: ii
integer :: j
integer :: ll
real(kind=PREC_DP) :: sum
!
II = 0 
!
DO I = 1, N 
   LL = INDX (I) 
   SUM = B (LL) 
   B (LL) = B (I) 
   IF (II.NE.0) THEN 
      DO J = II, I - 1 
         SUM = SUM - A (I, J) * B (J) 
      END DO 
   ELSEIF (SUM.NE.0.) THEN 
      II = I 
   ENDIF 
   B (I) = SUM 
END DO 
!
DO I = N, 1, - 1 
   SUM = B (I) 
   IF (I.LT.N) THEN 
      DO J = I + 1, N 
         SUM = SUM - A (I, J) * B (J) 
      END DO 
   ENDIF 
   B (I) = SUM / A (I, I) 
END DO 
!     RETURN 
!
END SUBROUTINE LUBKSB                         
!
!*******************************************************************************
!
SUBROUTINE LUDCMP (A, N, NP, INDX, D, ier) 
!
use precision_mod
!
implicit none
!
!IMPLICIT integer(i-n)
!IMPLICIT REAL (a-h, o-z)
!
integer                              , intent(in) :: N
integer                              , intent(in) :: NP
real(kind=PREC_DP), dimension(NP, NP), intent(inout) :: a
integer           , dimension(N)     , intent(out)   :: indx
real(kind=PREC_DP)                   , intent(out)   :: d
INTEGER                              , INTENT(out)   :: ier
!
real(KIND=PREC_DP), parameter :: TINY = 1.0D-20
!     PARAMETER (NMAX = 100, TINY = 1.0E-20) 
!      DIMENSION A (NP, NP), INDX (N), VV (NMAX) 
integer, parameter :: NMAX = 500
real(kind=PREC_DP), dimension(NMAX) :: vv
integer :: i
integer :: j
integer :: k
integer :: imax
real(kind=PREC_DP) :: aamax
real(kind=PREC_DP) :: dum
real(kind=PREC_DP) :: sum
!
imax = 0
D = 1. 
DO I = 1, N 
   AAMAX = 0. 
   DO J = 1, N 
      IF (ABS (A (I, J) ) .GT.AAMAX) AAMAX = ABS (A (I, J) ) 
   END DO 
   IF (AAMAX.EQ.0.) THEN !! PAUSE 'Singular matrix.' 
      ier = -61
      RETURN
   ENDIF
   VV (I) = 1. / AAMAX 
END DO 
!
DO J = 1, N 
   IF (J.GT.1) THEN 
      DO I = 1, J - 1 
         SUM = A (I, J) 
         IF (I.GT.1) THEN 
            DO K = 1, I - 1 
               SUM = SUM - A (I, K) * A (K, J) 
            END DO 
            A (I, J) = SUM 
         ENDIF 
      END DO 
   ENDIF 
   AAMAX = 0. 
   DO I = J, N 
      SUM = A (I, J) 
      IF (J.GT.1) THEN 
         DO K = 1, J - 1 
            SUM = SUM - A (I, K) * A (K, J) 
         END DO 
         A (I, J) = SUM 
      ENDIF 
      DUM = VV (I) * ABS (SUM) 
      IF (DUM.GE.AAMAX) THEN 
         IMAX = I 
         AAMAX = DUM 
      ENDIF 
   END DO 
   IF (J.NE.IMAX) THEN 
      DO K = 1, N 
         DUM = A (IMAX, K) 
         A (IMAX, K) = A (J, K) 
         A (J, K) = DUM 
      END DO 
      D = - D 
      VV (IMAX) = VV (J) 
   ENDIF 
   INDX (J) = IMAX 
   IF (J.NE.N) THEN 
      IF (A (J, J) .EQ.0.) A (J, J) = TINY 
      DUM = 1. / A (J, J) 
      DO I = J + 1, N 
         A (I, J) = A (I, J) * DUM 
      END DO 
   ENDIF 
END DO 
IF (A (N, N) .EQ.0.) A (N, N) = TINY 
!      RETURN 
!
END SUBROUTINE LUDCMP                         
!
!*******************************************************************************
!
end module kuplot_math_mod
