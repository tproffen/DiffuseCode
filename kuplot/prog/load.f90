!**************************************************************** 
!     Routines to allocate memory for data sets and load data sets      
!     using different file formats.                                     
!*****7**************************************************************** 
SUBROUTINE do_func (zeile, lp) 
!                                                                       
!     Create new data set from given function                           
!                                                                       
      USE ber_params_mod
      USE berechne_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE kuplot_fit6
USE precision_mod
USE str_comp_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 7) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw), cfkt, cdummy 
      INTEGER lp, lcfkt, lpara (maxw) 
      INTEGER ianz, i, ii, jj, kk, iref 
      INTEGER :: length
      INTEGER maxpkt, maxzz 
      REAL xstart, xend, xdelta 
      REAL ystart, yend, ydelta 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL df (maxpara) 
      REAL xx, yy, f, xkk 
      LOGICAL ffit
!                                                                       
!                                                                       
!------ space left for additional data set ?                            
!                                                                       
      IF (iz.gt.maxkurvtot) then 
         ier_num = - 1 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      lcfkt = 1
!                                                                       
      maxpkt = maxarray - offxy (iz - 1) 
      maxzz = maxarray - offz (iz - 1) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ianz.lt.1) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
      IF (str_comp (cpara (1) , 'fit', 3, lpara (1) , 3) ) then 
         ffit = .true. 
      ELSE 
         cfkt = cpara (1) (1:lpara (1) ) 
         lcfkt = lpara (1) 
         ffit = .false. 
      ENDIF 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
!                                                                       
!------ Use reference set                                               
!                                                                       
nianz: IF (ianz.eq.1) then 
   CALL ber_params (ianz, cpara, lpara, werte, maxw) 
   iref = nint (werte (1) ) 
!                                                                       
   IF (iref.le.0.or.iref.ge.iz) THEN 
      ier_num = - 4 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
!                                                                       
!  IF (lni (iref) ) THEN 
!     ier_num = - 4 
!     ier_typ = ER_APPL 
!     RETURN 
!  ENDIF 
!                                                                       
   ii = lenc(iref) 
!                                                                       
   IF (ii.le.0.or.ii.gt.maxpkt) THEN 
      ier_num = - 6 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
!                                                                       
!-------  KUPLOT fit function                                           
!                                                                       
   IF (ffit) THEN 
      IF(lni(iref)) THEN     ! 2d data set
         nx(iz) = nx(iref)
         ny(iz) = nx(iref)
         DO i = 1, nx (ikfit) * ny (ikfit)
            xx = REAL(i)
            CALL kupl_theory (xx, f, df, - i) 
            z (offz (iz - 1) + i) = f
         ENDDO
         DO ii = 1, nx (iz) 
            x (offxy (iz - 1) + ii) = x (offxy (iref - 1) + ii)
         ENDDO
         DO jj = 1, ny (iz) 
            y (offxy (iz - 1) + jj) = y (offxy (iref - 1) + jj)
         ENDDO
         fname (iz) = 'func.nipl' 
         fform (iz) = 'NI' 
         lni (iz) = .true. 
         lh5 (iz) = .FALSE.
         lenc(iz) = max (nx (iz), ny (iz) ) 
         offxy (iz) = offxy (iz - 1) + lenc(iz) 
         offz (iz) = offz (iz - 1) + nx (iz) * ny (iz) 
         iz = iz + 1 
      ELSE
         DO i = 1, ii 
            xx = x (offxy (iref - 1) + i) 
            CALL kupl_theory (xx, f, df, - i) 
            x (offxy (iz - 1) + i) = xx 
            y (offxy (iz - 1) + i) = f 
            dx (offxy (iz - 1) + i) = 0.0 
            dy (offxy (iz - 1) + i) = 0.0 
         ENDDO 
         fname (iz) = 'func.xy' 
         fform (iz) = 'XY' 
         lni (iz) = .false. 
         lh5 (iz) = .FALSE.
         lenc(iz) = ii 
         offxy (iz) = offxy (iz - 1) + lenc(iz) 
         offz (iz) = offz (iz - 1) 
         iz = iz + 1 
      ENDIF
!                                                                       
!------- User defined function                                          
!                                                                       
   ELSE 
      IF(lni(iref)) THEN     ! 2d data set
         nx(iz) = nx(iref)
         ny(iz) = nx(iref)
         DO ii = 1, nx (iz) 
            rpara (0) = x (offxy (iref - 1) + ii) 
            DO jj = 1, ny (iz) 
               cdummy = '('//cfkt (1:lcfkt) //')' 
               length = lcfkt + 2
               rpara (1) = y (offxy (iref - 1) + jj) 
               z (offz (iz - 1) + (ii - 1) * ny (iz) + jj) =  &
                 berechne(cdummy, length)
               x (offxy (iz - 1) + ii) = rpara (0) 
               y (offxy (iz - 1) + jj) = rpara (1) 
               dx (offxy (iz - 1) + ii) = 0.0 
               dy (offxy (iz - 1) + jj) = 0.0 
            ENDDO 
         ENDDO 
!
         fname (iz) = 'func.nipl' 
         fform (iz) = 'NI' 
         lni (iz) = .true. 
         lh5 (iz) = .FALSE.
         lenc(iz) = max (nx (iz), ny (iz) ) 
         offxy (iz) = offxy (iz - 1) + lenc(iz) 
         offz (iz) = offz (iz - 1) + nx (iz) * ny (iz) 
         iz = iz + 1 
      ELSE
         DO i = 1, ii 
            rpara (0) = x (offxy (iref - 1) + i) 
            cdummy = '('//cfkt (1:lcfkt) //')' 
            length = lcfkt + 2
            x (offxy (iz - 1) + i) = rpara (0) 
            y (offxy (iz - 1) + i) = berechne (cdummy, length)
            dx (offxy (iz - 1) + i) = 0.0 
            dy (offxy (iz - 1) + i) = 0.0 
         ENDDO 
         fname (iz) = 'func.xy' 
         fform (iz) = 'XY' 
         lni (iz) = .false. 
         lh5 (iz) = .FALSE.
         lenc(iz) = ii 
         offxy (iz) = offxy (iz - 1) + lenc(iz) 
         offz (iz) = offz (iz - 1) 
         iz = iz + 1 
      ENDIF 
   ENDIF 
   CALL show_data (iz - 1) 
!                                                                       
!------ 1D data set                                                     
!                                                                       
ELSEIF (ianz.eq.3) THEN  nianz
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         xstart = werte (1) 
         xend = werte (2) 
         xdelta = werte (3) 
         ii = nint ( (xend-xstart) / xdelta) + 1 
         IF (ii.gt.0.and.ii.le.maxpkt) then 
!                                                                       
!------- -- KUPLOT fit function                                         
!                                                                       
            IF (ffit) then 
               DO i = 1, ii 
               xx = xstart + (i - 1) * xdelta 
               CALL kupl_theory (xx, f, df, - i) 
               x (offxy (iz - 1) + i) = xx 
               y (offxy (iz - 1) + i) = f 
               dx (offxy (iz - 1) + i) = 0.0 
               dy (offxy (iz - 1) + i) = 0.0 
               ENDDO 
!                                                                       
!------- -- User defined function                                       
!                                                                       
            ELSE 
               DO i = 1, ii 
               rpara (0) = xstart + (i - 1) * xdelta 
               cdummy = '('//cfkt (1:lcfkt) //')' 
               length = lcfkt + 2
               x (offxy (iz - 1) + i) = rpara (0) 
               y (offxy (iz - 1) + i) = berechne (cdummy, length)
               dx (offxy (iz - 1) + i) = 0.0 
               dy (offxy (iz - 1) + i) = 0.0 
               ENDDO 
            ENDIF 
            fname (iz) = 'func.xy' 
            fform (iz) = 'XY' 
            lni (iz) = .false. 
            lh5 (iz) = .FALSE.
            lenc(iz) = ii 
            offxy (iz) = offxy (iz - 1) + lenc(iz) 
            offz (iz) = offz (iz - 1) 
            iz = iz + 1 
            CALL show_data (iz - 1) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ 2D data set                                                     
!                                                                       
ELSEIF (ianz.eq.6) THEN nianz
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         xstart = werte (1) 
         xend = werte (2) 
         xdelta = werte (3) 
         ystart = werte (4) 
         yend = werte (5) 
         ydelta = werte (6) 
         nx (iz) = nint ( (xend-xstart) / xdelta) + 1 
         ny (iz) = nint ( (yend-ystart) / ydelta) + 1 
         IF (nx (iz) .gt.0.and.ny (iz) .gt.0.and.nx (iz) * ny (iz)      &
         .le.maxzz.and.max (nx (iz), ny (iz) ) .le.maxpkt) then         
!                                                                       
!------- -- KUPLOT fit function                                         
!                                                                       
            IF (ffit) then 
               DO ii = 1, nx (iz) 
               xx = xstart + (ii - 1) * xdelta 
               DO jj = 1, ny (iz) 
               kk = (ii - 1) * ny (iz) + jj 
               xkk = REAL(kk) 
               yy = ystart + (jj - 1) * ydelta 
               CALL kupl_theory (xkk, f, df, - kk) 
               z (offz (iz - 1) + (ii - 1) * ny (iz) + jj) = f 
               x (offxy (iz - 1) + ii) = xx 
               y (offxy (iz - 1) + jj) = yy 
               dx (offxy (iz - 1) + ii) = 0.0 
               dy (offxy (iz - 1) + jj) = 0.0 
               ENDDO 
               ENDDO 
!                                                                       
!------- -- User defined function                                       
!                                                                       
            ELSE 
               DO ii = 1, nx (iz) 
               rpara (0) = xstart + (ii - 1) * xdelta 
               DO jj = 1, ny (iz) 
               cdummy = '('//cfkt (1:lcfkt) //')' 
               length = lcfkt + 2
               rpara (1) = ystart + (jj - 1) * ydelta 
               z (offz (iz - 1) + (ii - 1) * ny (iz) + jj) =  &
                 berechne(cdummy, length)
               x (offxy (iz - 1) + ii) = rpara (0) 
               y (offxy (iz - 1) + jj) = rpara (1) 
               dx (offxy (iz - 1) + ii) = 0.0 
               dy (offxy (iz - 1) + jj) = 0.0 
               ENDDO 
               ENDDO 
            ENDIF 
!                                                                       
            fname (iz) = 'func.nipl' 
            fform (iz) = 'NI' 
            lni (iz) = .true. 
            lh5 (iz) = .FALSE.
            lenc(iz) = max (nx (iz), ny (iz) ) 
            offxy (iz) = offxy (iz - 1) + lenc(iz) 
            offz (iz) = offz (iz - 1) + nx (iz) * ny (iz) 
            iz = iz + 1 
            CALL show_data (iz - 1) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
ELSE  nianz
         ier_num = - 6 
         ier_typ = ER_COMM 
ENDIF  nianz
!                                                                       
      END SUBROUTINE do_func                        
!*****7**************************************************************** 
      SUBROUTINE do_allocate (zeile, lp, lecho) 
!                                                                       
!     Allocate space for new data set                                   
!                                                                       
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE lib_errlist_func
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: zeile 
INTEGER         , INTENT(INOUT) :: lp
LOGICAL         , INTENT(IN)    :: lecho
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(40) cdummy 
      INTEGER lpara (maxw) 
      INTEGER ianz, i, ii, jj 
      INTEGER maxpkt, maxzz 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL no_error 
!                                                                       
!------ space left for additional data set ?                            
!                                                                       
      IF (iz.gt.maxkurvtot) then 
         ier_num = - 1 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      maxpkt = maxarray - offxy (iz - 1) 
      maxzz = maxarray - offz (iz - 1) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.ne.2.and.ianz.ne.3) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
                                                                        
      cdummy = cpara (1) (1:MIN(40,LEN_TRIM(cpara(1))))
      CALL del_params (1, ianz, cpara, lpara, maxw) 
!                                                                       
!------ 1D data set                                                     
!                                                                       
      IF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         ii = nint (werte (1) ) 
         IF (ii.gt.0.and.ii.le.maxpkt) then 
            DO i = 1, ii 
            x (offxy (iz - 1) + i) = 0.0 
            y (offxy (iz - 1) + i) = 0.0 
            dx (offxy (iz - 1) + i) = 0.0 
            dy (offxy (iz - 1) + i) = 0.0 
            ENDDO 
            fname (iz) = cdummy 
            fform (iz) = 'XY' 
            lni (iz) = .false. 
            lh5 (iz) = .FALSE.
            lenc(iz) = ii 
            offxy (iz) = offxy (iz - 1) + lenc(iz) 
            offz (iz) = offz (iz - 1) 
            iz = iz + 1 
            IF(lecho) CALL show_data (iz - 1) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ 2D data set                                                     
!                                                                       
      ELSEIF (ianz.eq.2) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         nx (iz) = nint (werte (1) ) 
         ny (iz) = nint (werte (2) ) 
         IF (nx (iz) .gt.0.and.ny (iz) .gt.0.and.nx (iz) * ny (iz)      &
         .le.maxzz.and.max (nx (iz), ny (iz) ) .le.maxpkt) then         
            DO ii = 1, nx (iz) 
            DO jj = 1, ny (iz) 
            z (offz (iz - 1) + (ii - 1) * ny (iz) + jj) = 0.0 
            x (offxy (iz - 1) + ii) = 0.0 
            y (offxy (iz - 1) + jj) = 0.0 
            dx (offxy (iz - 1) + ii) = 0.0 
            dy (offxy (iz - 1) + jj) = 0.0 
            ENDDO 
            ENDDO 
            fname (iz) = cdummy 
            fform (iz) = 'NI' 
            lni (iz) = .true. 
            lh5 (iz) = .FALSE.
            lenc(iz) = max (nx (iz), ny (iz) ) 
            offxy (iz) = offxy (iz - 1) + lenc(iz) 
            offz (iz) = offz (iz - 1) + nx (iz) * ny (iz) 
            iz = iz + 1 
            IF(lecho) CALL show_data (iz - 1) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_allocate                    
!*****7**************************************************************** 
      SUBROUTINE do_load (string, laenge, lecho) 
!                                                                       
!     Load various file formats                                         
!                                                                       
      USE ber_params_mod
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE times_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE kuplot_load_h5
USE lib_errlist_func
USE lib_length
USE precision_mod
      USE take_param_mod
USE str_comp_mod
      USE string_convert_mod
USE support_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw, ifil, iwgb 
      PARAMETER (maxw = 100) 
      PARAMETER (ifil = 44) 
      PARAMETER (iwgb = 45) 
!                                                                       
      CHARACTER ( * ) string 
!
      INTEGER, INTENT(INOUT) :: laenge
      LOGICAL, INTENT(IN) :: lecho          ! Show extrema after load
!
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(LEN=PREC_STRING) :: wname, cdummy 
      CHARACTER(4) unter 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ii, ll, ianz, istr, nfile 
      LOGICAL zz_mod 
!                                                                       
INTEGER ifiles
!
INTEGER, PARAMETER :: NOPTIONAL = 7
!INTEGER, PARAMETER :: O_SKIP      = 1
!INTEGER, PARAMETER :: O_COLX      = 2
!INTEGER, PARAMETER :: O_COLY      = 3
!INTEGER, PARAMETER :: O_COLDX     = 4
!INTEGER, PARAMETER :: O_COLDY     = 5
INTEGER, PARAMETER :: O_LAYER     = 6
!INTEGER, PARAMETER :: O_SEPARATOR = 7
CHARACTER(LEN=          9), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 5 ! Number of values to calculate 
!
DATA oname  / 'skip', 'colx',  'coly',  'coldx', 'coldy', 'layer', 'separator'  /
DATA loname /  4    ,  4    ,   4    ,   5     ,  5     ,  5     ,  9           /
opara  =  (/ '25.000', '1.0000', '2.0000', '0.0000', '0.0000', 'middle', ';     ' /)   ! Always provide fresh default values
lopara =  (/  6,        6,        6      ,  6      ,  6      ,  6      ,  6       /)
owerte =  (/ 25.0,      1.0,      2.0    ,  0.0    ,  0.0    ,  1.0    ,  0.0     /)
!
!
!                                                                       
      istr = 1 
      CALL no_error 
      CALL get_params (string, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.ne.0) return 
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      IF(ier_num/=0) RETURN
      IF (ianz.ge.2) then 
!                                                                       
!------- -get file typ                                                  
!                                                                       
         unter = cpara (1) (1:4) 
         CALL do_cap (unter) 
!                                                                       
!------- -build first filename                                          
!                                                                       
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 2) 
         IF (ier_num.ne.0) return 
!                                                                       
         nfile = ifiles (cpara (2), lpara (2) ) 
         IF (nfile.eq.0) then 
            ier_num = - 1 
            ier_typ = ER_IO 
            RETURN 
         ENDIF 
!
         IF(unter=='H5') THEN    ! HDF5 file
            CALL hdf5_read(cpara(2),lpara(2), O_LAYER, NOPTIONAL, opara,        &
                           lopara, lpresent, owerte ,                           &
                           MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny,    &
                           xmin, xmax, ymin, ymax, offxy, offz, lni, lh5, lenc, &
                           ier_num, ier_typ, UBOUND(ier_msg,1), ier_msg, ER_APPL, &
                           ER_IO, output_io)
            RETURN
         ENDIF
!                                                                       
         IF (ianz.eq.3.and.unter.ne.'SP') then 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 3) 
            IF (ier_num.ne.0) return 
            istr = 2 
            wname = cpara (3) 
         ELSE 
            istr = 1 
            wname = ' ' 
         ENDIF 
!                                                                       
!------- -for format 'ZZ' or 'DE' get grid size                         
!                                                                       
         IF ( (unter.eq.'ZZ'.or.unter.eq.'DE'.or.unter.eq.'MP') ) then 
            IF (ianz.eq.4.or.ianz.eq.5.or.ianz.eq.8) then 
               zz_mod = .false. 
               IF (ianz.eq.5) then 
                  IF (str_comp (cpara (ianz) , 'mod', 3, lpara (ianz) , &
                  3) ) then                                             
                     zz_mod = .true. 
                     ianz = 4 
                  ENDIF 
               ENDIF 
               CALL del_params (2, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
!                                                                       
!------ --Other additional parameters                                   
!                                                                       
         ELSEIF (unter.eq.'MC'.or.unter.eq.'XY'.or.unter.eq.'4D') then 
            IF (ianz.gt.2) then 
               CALL del_params (2, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
            ELSE 
               ianz = 0 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!
!                                                                       
!------ Here we open all data files                                     
!                                                                       
      DO ii = 1, nfile 
      IF (iz.gt.maxkurvtot) then 
         ier_num = - 1 
         ier_typ = ER_APPL 
         CALL freefiles 
         RETURN 
      ENDIF 
!                                                                       
      cdummy = '' 
      ll = 200 
      CALL getfile (cdummy, ll, ii) 
      fname (iz) = cdummy (1:ll)
      IF(lecho) WRITE (output_io, 2000) fname (iz) (1:len_str (fname (iz) ) ) 
!                                                                       
      IF (unter.ne.'MP') then 
         CALL oeffne (ifil, fname (iz) , 'old') 
      ELSE 
         OPEN (ifil, file = fname (iz) , status = 'unknown', access =   &
         'direct', form = 'unformatted', recl = 16)                     
         CALL file_info (ifil) 
      ENDIF 
      IF (ier_num.ne.0) return 
      IF (                                                              &
      istr.eq.2.and.unter.ne.'ZZ'.and.unter.ne.'DE'.and.unter.ne.'MP')  &
      then                                                              
         CALL oeffne (iwgb, wname, 'old') 
         IF (ier_num.ne.0) then 
            istr = 1 
            ier_num = 0 
            ier_typ = ER_NONE 
         ENDIF 
      ENDIF 
!                                                                       
!------ - Has the file a prepended History ?                            
!                                                                       
      IF (unter.ne.'MP') then 
         CALL no_error 
         CALL extract_hist (ifil) 
         IF (ier_num.ne.0) return 
      ENDIF 
!                                                                       
!------ - Check file format                                             
!                                                                       
      fform (iz) = unter (1:2)
      IF (unter.eq.'XY') then 
         CALL read_xy (ifil, ianz, werte, maxw) 
      ELSEIF (unter.eq.'XX') then 
         CALL read_xx (ifil) 
      ELSEIF (unter.eq.'CR') then 
         CALL read_crystal (ifil) 
      ELSEIF (unter.eq.'DE') then 
         CALL read_xyz (ifil, ianz, werte, maxw, .true., zz_mod) 
      ELSEIF (unter.eq.'ZZ') then 
         CALL read_xyz (ifil, ianz, werte, maxw, .false., zz_mod) 
      ELSEIF (unter.eq.'4D') then 
         CALL read_4d (ifil, ianz, werte, maxw) 
      ELSEIF (unter.eq.'MC') then 
         CALL read_mca_single (ifil, ianz, werte, maxw) 
      ELSEIF (unter.eq.'MP') then 
         CALL read_ext (ifil, ianz, werte, maxw, .false., zz_mod) 
      ELSEIF (unter.eq.'Z1') then 
         CALL read_z (ifil, 1) 
      ELSEIF (unter.eq.'Z2') then 
         CALL read_z (ifil, 2) 
      ELSEIF (unter.eq.'MA') then 
         CALL read_marker (ifil) 
      ELSEIF (unter.eq.'NI') then 
         CALL read_nipl (ifil, iwgb, istr, .false.) 
      ELSEIF (unter.eq.'PG') then 
         CALL read_nipl (ifil, iwgb, istr, .true.) 
      ELSEIF (unter.eq.'T1') then 
         CALL read_mess (ifil, 4, 8, '2-theta [deg.]', 'counts z1') 
      ELSEIF (unter.eq.'T2') then 
         CALL read_mess (ifil, 4, 9, '2-theta [deg.]', 'counts z2') 
      ELSEIF (unter.eq.'C1') then 
         CALL read_mess (ifil, 5, 8, 'chi [deg.]', 'counts z1') 
      ELSEIF (unter.eq.'C2') then 
         CALL read_mess (ifil, 5, 9, 'chi [deg.]', 'counts z2') 
      ELSEIF (unter.eq.'P1') then 
         CALL read_mess (ifil, 6, 8, 'phi [deg.]', 'counts z1') 
      ELSEIF (unter.eq.'P2') then 
         CALL read_mess (ifil, 6, 9, 'phi [deg.]', 'counts z2') 
      ELSEIF (unter.eq.'O1') then 
         CALL read_mess (ifil, 7, 8, 'omega [deg.]', 'counts z1') 
      ELSEIF (unter.eq.'O2') then 
         CALL read_mess (ifil, 7, 9, 'omega [deg.]', 'counts z2') 
      ELSEIF (unter.eq.'H1') then 
         CALL read_mess (ifil, 1, 8, 'h [r.l.u.]', 'counts z1') 
      ELSEIF (unter.eq.'H2') then 
         CALL read_mess (ifil, 1, 9, 'h [r.l.u.]', 'counts z2') 
      ELSEIF (unter.eq.'K1') then 
         CALL read_mess (ifil, 2, 8, 'k [r.l.u.]', 'counts z1') 
      ELSEIF (unter.eq.'K2') then 
         CALL read_mess (ifil, 2, 9, 'k [r.l.u.]', 'counts z2') 
      ELSEIF (unter.eq.'L1') then 
         CALL read_mess (ifil, 3, 8, 'l [r.l.u.]', 'counts z1') 
      ELSEIF (unter.eq.'L2') then 
         CALL read_mess (ifil, 3, 9, 'l [r.l.u.]', 'counts z2') 
      ELSEIF (unter.eq.'TE') then 
         CALL read_mess (ifil, 0, 11, 'data point', 'temperature [k]') 
      ELSEIF (unter.eq.'DM') then 
         CALL read_mess (ifil, 0, 12, 'data point', 'value DMM') 
      ELSEIF (unter.eq.'ZE') then 
         CALL read_mess (ifil, 0, 10, 'data point', 'time [sec]') 
      ELSEIF (unter.eq.'SC'.or.unter.eq.'ST'.or.unter.eq.'SMCA') then 
         IF (unter.eq.'SC') then 
            CALL read_spec (ifil, ianz, cpara, lpara, werte, maxw, 0) 
         ELSEIF (unter.eq.'ST') then 
            CALL read_spec (ifil, ianz, cpara, lpara, werte, maxw, 1) 
         ELSEIF (unter.eq.'SMCA') then 
            CALL read_spec (ifil, ianz, cpara, lpara, werte, maxw, 2) 
         ENDIF 
      ELSEIF (unter.eq.'GS') then 
         CALL read_gsas (ifil, ianz, cpara, lpara, werte, maxw) 
      ELSEIF (unter.eq.'SM') then 
         CALL read_simref (ifil) 
      ELSEIF (unter.eq.'CSV') then 
         CALL do_read_csv(maxw, ianz, cpara, lpara ,   &
              NOPTIONAL, opara, lopara, owerte, ncalc , ifil, lecho)
      ELSE 
         ier_num = - 2 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      CLOSE (ifil) 
      IF (istr.eq.2) close (iwgb) 
      ENDDO 
      CALL freefiles 
      RETURN 
!                                                                       
      CLOSE (ifil) 
      IF (istr.eq.2) close (iwgb) 
!                                                                       
      ier_num = - 40 
      ier_typ = ER_APPL 
!                                                                       
 2000 FORMAT    (1x,'Reading file ',a,' ..') 
!                                                                       
      END SUBROUTINE do_load                        
!*****7**************************************************************** 
      SUBROUTINE extract_hist (ifil) 
!                                                                       
!     Extract history information if present                            
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ifil 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: line 
      INTEGER is
!                                                                       
!                                                                       
      READ (ifil, 5000, end = 40, err = 44) line 
      IF (line (1:7) .eq.'History') then 
         res_para (0) = 0 
         WRITE (output_io, 1000) 
   22    CONTINUE 
         READ (ifil, 5000) line 
!                                                                       
         is = index (line, 'runTitle=') 
         IF (is.ne.0) then 
            titel (iwin, iframe, 1) = line (is + 9:len_str (line) ) 
         ENDIF 
!                                                                       
         CALL extract_key (line, 'temp=') 
         CALL extract_key (line, 'Qmax=') 
         IF (line (1:16) .eq.'##### start data') goto 33 
         GOTO 22 
      ELSE 
         BACKSPACE (ifil) 
      ENDIF 
!                                                                       
   33 CONTINUE 
      RETURN 
   40 CONTINUE 
      ier_num = - 6 
      ier_typ = ER_IO 
      RETURN 
   44 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
 1000 FORMAT    ( ' History information found ...') 
 5000 FORMAT    (a) 
!                                                                       
      END SUBROUTINE extract_hist                   
!*****7*****************************************************************
      SUBROUTINE extract_key (line, key) 
!+                                                                      
!     Gets numbers from history part of data file                       
!-                                                                      
      USE param_mod 
      USE kuplot_config 
USE lib_length
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) line, key 
      INTEGER is, ie, ll, lk 
!                                                                       
!                                                                       
      ll = len_str (line) 
      lk = len_str (key) 
!                                                                       
      is = index (line, key (1:lk) ) 
      IF (is.ne.0) then 
         is = is + lk 
         DO while (line (is:is) .eq.' ') 
         is = is + 1 
         ENDDO 
         ie = index (line (is:ll) , ' ') 
         IF (ie.eq.0) ie = ll 
         res_para (0) = res_para (0) + 1 
         READ (line (is:is + ie), * ) res_para (NINT(res_para (0)) ) 
      ENDIF 
!                                                                       
      END SUBROUTINE extract_key                    
!*****7**************************************************************** 
      SUBROUTINE read_ext (ifil, ianz, werte, maxw, dens, zz_mod) 
!+                                                                      
!     Load MPAUS extract files ..                                       
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
!                                                                       
      REAL(KIND=PREC_DP) :: werte (maxw) 
!DBG      real xw,yw,zw                                                 
      REAL xw, yw, ow 
      INTEGER zw 
      REAL deltax, deltay, dxx, dyy 
      INTEGER irec 
      INTEGER np (maxarray) 
      INTEGER ifil, ixx, iyy, izeig, i 
      INTEGER maxpkt, maxzz, ianz 
      LOGICAL dens, zz_mod 
!                                                                       
!------ get parameters                                                  
!                                                                       
      fform (iz) = 'NI' 
!                                                                       
      maxpkt = maxarray - offxy (iz - 1) 
      maxzz = maxarray - offz (iz - 1) 
      deltax = werte (1) 
      deltay = werte (2) 
!                                                                       
!     get number of records                                             
!                                                                       
      READ (ifil, rec = 1) irec 
!                                                                       
!------ get dimensions of file                                          
!                                                                       
      IF (ianz.eq.6) then 
         xmin (iz) = werte (3) 
         xmax (iz) = werte (4) 
         ymin (iz) = werte (5) 
         ymax (iz) = werte (6) 
      ELSE 
         READ (ifil, rec = 2) xw, yw, zw, ow 
         xmin (iz) = xw 
         xmax (iz) = xw 
         ymin (iz) = yw 
         ymax (iz) = yw 
         DO i = 3, irec 
         READ (ifil, rec = i) xw, yw, zw, ow 
         xmin (iz) = min (xmin (iz), xw) 
         xmax (iz) = max (xmax (iz), xw) 
         ymin (iz) = min (ymin (iz), yw) 
         ymax (iz) = max (ymax (iz), yw) 
         ENDDO 
      ENDIF 
!                                                                       
!     optionally round xmin... to integer multiples of delta            
!                                                                       
      IF (zz_mod) then 
         xmin (iz) = REAL(nint (xmin (iz) / deltax) ) * deltax 
         xmax (iz) = REAL(nint (xmax (iz) / deltax) ) * deltax 
         ymin (iz) = REAL(nint (ymin (iz) / deltay) ) * deltay 
         ymax (iz) = REAL(nint (ymax (iz) / deltay) ) * deltay 
      ENDIF 
!                                                                       
!------ check dimensions                                                
!                                                                       
      IF (xmin (iz) .eq.xmax (iz) .or.ymin (iz) .eq.ymax (iz) ) then 
         ier_num = - 29 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ get grid size                                                   
!                                                                       
      nx (iz) = nint ( (xmax (iz) - xmin (iz) ) / deltax) + 1 
      ny (iz) = nint ( (ymax (iz) - ymin (iz) ) / deltay) + 1 
!                                                                       
!-------check array size                                                
!                                                                       
      IF ( (nx (iz) * ny (iz) .gt.maxzz) .or. (max (nx (iz), ny (iz) )  &
      .gt.maxpkt) ) then                                                
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ zero some arrays                                                
!                                                                       
      DO ixx = 1, nx (iz) 
      DO iyy = 1, ny (iz) 
      izeig = offz (iz - 1) + (ixx - 1) * ny (iz) + iyy 
      z (izeig) = 0.0 
      np (izeig) = 0
      ENDDO 
      ENDDO 
!                                                                       
!------ read file and bin on specified grid                             
!                                                                       
      dxx = deltax 
      dyy = deltay 
!                                                                       
      DO i = 2, irec 
      READ (ifil, rec = i) xw, yw, zw 
      IF (zw.ne. - 9999.0) then 
         ixx = nint ( (xw - xmin (iz) ) / dxx) + 1 
         iyy = nint ( (yw - ymin (iz) ) / dyy) + 1 
         IF (ixx.gt.0.and.ixx.le.nx (iz) .and.iyy.gt.0.and.iyy.le.ny (  &
         iz) ) then                                                     
            izeig = offz (iz - 1) + (ixx - 1) * ny (iz) + iyy 
            z (izeig) = z (izeig) + REAL(zw) 
            np (izeig) = np (izeig) + 1 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!------ normalise Z                                                     
!                                                                       
      DO ixx = 1, nx (iz) 
      DO iyy = 1, ny (iz) 
      izeig = offz (iz - 1) + (ixx - 1) * ny (iz) + iyy 
      IF (.not.dens) then 
         IF (np (izeig) .eq.0) then 
            z (izeig) = - 9999.0 
         ELSE 
            z (izeig) = z (izeig) / REAL(np (izeig) ) 
         ENDIF 
      ENDIF 
      zmax (iz) = max (z (izeig), zmax (iz) ) 
      zmin (iz) = min (z (izeig), zmin (iz) ) 
      ENDDO 
      ENDDO 
!                                                                       
!------ fill X,Y arrays                                                 
!                                                                       
      DO i = 1, nx (iz) 
      x (offxy (iz - 1) + i) = xmin (iz) + (i - 1) * dxx 
      ENDDO 
      DO i = 1, ny (iz) 
      y (offxy (iz - 1) + i) = ymin (iz) + (i - 1) * dyy 
      ENDDO 
!                                                                       
!------ set other parameters                                            
!                                                                       
      lni (iz) = .true. 
      lh5 (iz) = .FALSE.
      lenc(iz) = max (nx (iz), ny (iz) ) 
      offxy (iz) = offxy (iz - 1) + lenc(iz) 
      offz (iz) = offz (iz - 1) + nx (iz) * ny (iz) 
      iz = iz + 1 
      CALL show_data (iz - 1) 
!                                                                       
      END SUBROUTINE read_ext                       
!*****7**************************************************************** 
      SUBROUTINE read_xyz (ifil, ianz, werte, maxw, dens, zz_mod) 
!+                                                                      
!     Load xyz files ..                                                 
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
!                                                                       
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL xw, yw, zw 
      REAL deltax, deltay, dxx, dyy 
      INTEGER np (maxarray) 
      INTEGER ifil, ixx, iyy, izeig, i 
      INTEGER maxpkt, maxzz, ianz 
      LOGICAL dens, zz_mod 
!                                                                       
!------ get parameters                                                  
!                                                                       
      fform (iz) = 'NI' 
!                                                                       
      maxpkt = maxarray - offxy (iz - 1) 
      maxzz = maxarray - offz (iz - 1) 
      deltax = werte (1) 
      deltay = werte (2) 
!                                                                       
!------ get dimensions of file                                          
!                                                                       
      IF (ianz.eq.6) then 
         xmin (iz) = werte (3) 
         xmax (iz) = werte (4) 
         ymin (iz) = werte (5) 
         ymax (iz) = werte (6) 
      ELSE 
         READ (ifil, *, end = 20) xw, yw, zw 
         xmin (iz) = xw 
         xmax (iz) = xw 
         ymin (iz) = yw 
         ymax (iz) = yw 
   10    CONTINUE 
         READ (ifil, *, end = 20) xw, yw, zw 
         xmin (iz) = min (xmin (iz), xw) 
         xmax (iz) = max (xmax (iz), xw) 
         ymin (iz) = min (ymin (iz), yw) 
         ymax (iz) = max (ymax (iz), yw) 
         GOTO 10 
   20    CONTINUE 
         REWIND (ifil) 
      ENDIF 
!                                                                       
!     optionally round xmin... to integer multiples of delta            
!                                                                       
      IF (zz_mod) then 
         xmin (iz) = REAL(nint (xmin (iz) / deltax) ) * deltax 
         xmax (iz) = REAL(nint (xmax (iz) / deltax) ) * deltax 
         ymin (iz) = REAL(nint (ymin (iz) / deltay) ) * deltay 
         ymax (iz) = REAL(nint (ymax (iz) / deltay) ) * deltay 
      ENDIF 
!                                                                       
!------ check dimensions                                                
!                                                                       
      IF (xmin (iz) .eq.xmax (iz) .or.ymin (iz) .eq.ymax (iz) ) then 
         ier_num = - 29 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ get grid size                                                   
!                                                                       
      nx (iz) = nint ( (xmax (iz) - xmin (iz) ) / deltax) + 1 
      ny (iz) = nint ( (ymax (iz) - ymin (iz) ) / deltay) + 1 
!                                                                       
!-------check array size                                                
!                                                                       
      IF ( (nx (iz) * ny (iz) .gt.maxzz) .or. (max (nx (iz), ny (iz) )  &
      .gt.maxpkt) ) then                                                
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ zero some arrays                                                
!                                                                       
      DO ixx = 1, nx (iz) 
      DO iyy = 1, ny (iz) 
      izeig = offz (iz - 1) + (ixx - 1) * ny (iz) + iyy 
      z (izeig) = 0.0 
      np (izeig) = 0
      ENDDO 
      ENDDO 
!                                                                       
!------ read file and bin on specified grid                             
!                                                                       
      dxx = deltax 
      dyy = deltay 
   30 CONTINUE 
      IF (dens) then 
         READ (ifil, *, end = 40) xw, yw 
         zw = 1.0 
      ELSE 
         READ (ifil, *, end = 40) xw, yw, zw 
      ENDIF 
      IF (zw.ne. - 9999.0) then 
         ixx = nint ( (xw - xmin (iz) ) / dxx) + 1 
         iyy = nint ( (yw - ymin (iz) ) / dyy) + 1 
         IF (ixx.gt.0.and.ixx.le.nx (iz) .and.iyy.gt.0.and.iyy.le.ny (  &
         iz) ) then                                                     
            izeig = offz (iz - 1) + (ixx - 1) * ny (iz) + iyy 
            z (izeig) = z (izeig) + zw 
            np (izeig) = np (izeig) + 1 
         ENDIF 
      ENDIF 
      GOTO 30 
   40 CONTINUE 
!                                                                       
!------ normalise Z                                                     
!                                                                       
      DO ixx = 1, nx (iz) 
      DO iyy = 1, ny (iz) 
      izeig = offz (iz - 1) + (ixx - 1) * ny (iz) + iyy 
      IF (.not.dens) then 
         IF (np (izeig) .eq.0) then 
            z (izeig) = - 9999.0 
         ELSE 
            z (izeig) = z (izeig) / REAL(np (izeig) ) 
         ENDIF 
      ENDIF 
      zmax (iz) = max (z (izeig), zmax (iz) ) 
      zmin (iz) = min (z (izeig), zmin (iz) ) 
      ENDDO 
      ENDDO 
!                                                                       
!------ fill X,Y arrays                                                 
!                                                                       
      DO i = 1, nx (iz) 
      x (offxy (iz - 1) + i) = xmin (iz) + (i - 1) * dxx 
      ENDDO 
      DO i = 1, ny (iz) 
      y (offxy (iz - 1) + i) = ymin (iz) + (i - 1) * dyy 
      ENDDO 
!                                                                       
!------ set other parameters                                            
!                                                                       
      lni (iz) = .true. 
      lh5 (iz) = .FALSE.
      lenc(iz) = max (nx (iz), ny (iz) ) 
      offxy (iz) = offxy (iz - 1) + lenc(iz) 
      offz (iz) = offz (iz - 1) + nx (iz) * ny (iz) 
      iz = iz + 1 
      CALL show_data (iz - 1) 
!                                                                       
      END SUBROUTINE read_xyz                       
!*****7**************************************************************** 
      SUBROUTINE read_4d (ifil, ianz, werte, maxw) 
!+                                                                      
!     Load cut from xyzv (4D) files                                     
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
!                                                                       
      REAL val (maxz) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL range (2, 3) 
      REAL dxx, dyy 
      INTEGER np (3), icut (3), iii (3) 
      INTEGER maxpkt, maxzz, ianz 
      INTEGER ifil, ival 
      INTEGER i, j, izeig 
      INTEGER ixx, iyy, izz, jxx, jyy 
!                                                                       
      WRITE (output_io, 1000) 
 1000 FORMAT    (1x,'Be patient, this is really SLOW ..') 
!                                                                       
!------ get parameters                                                  
!                                                                       
      fform (iz) = 'NI' 
!                                                                       
      maxpkt = maxarray - offxy (iz - 1) 
      maxzz = maxarray - offz (iz - 1) 
      IF (ianz.ne.4) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
!------ get header information                                          
!                                                                       
      READ (ifil, * ) np 
      READ (ifil, * ) ( (range (i, j), i = 1, 2), j = 1, 3) 
      DO i = 1, 3 
      icut (i) = nint (werte (i) ) 
      ENDDO 
      ival = nint (werte (4) ) 
!                                                                       
      IF (icut (1) .lt.1.or.icut (1) .gt.3.or.icut (2) .lt.1.or.icut (2)&
      .gt.3.or.icut (3) .lt.1.or.icut (3)                               &
      .gt.3.or.ival.lt.1.or.ival.gt.np (icut (3) ) ) then               
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      nx (iz) = np (icut (1) ) 
      ny (iz) = np (icut (2) ) 
      xmin (iz) = range (1, icut (1) ) 
      xmax (iz) = range (2, icut (1) ) 
      ymin (iz) = range (1, icut (2) ) 
      ymax (iz) = range (2, icut (2) ) 
!                                                                       
      dxx = (xmax (iz) - xmin (iz) ) / REAL(nx (iz) - 1) 
      dyy = (ymax (iz) - ymin (iz) ) / REAL(ny (iz) - 1) 
!                                                                       
!-------check array size                                                
!                                                                       
      IF ( (nx (iz) * ny (iz) .gt.maxzz) .or. (max (nx (iz), ny (iz) )  &
      .gt.maxpkt) ) then                                                
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ read the data now                                               
!                                                                       
      DO ixx = 1, np (1) 
      DO iyy = 1, np (2) 
      iii (1) = ixx 
      iii (2) = iyy 
      READ (ifil, * ) (val (izz), izz = 1, np (3) ) 
      DO izz = 1, np (3) 
      IF ( (icut (3) .eq.1.and.ival.eq.ixx) .or. (icut (3)              &
      .eq.2.and.ival.eq.iyy) .or. (icut (3) .eq.3.and.ival.eq.izz) )    &
      then                                                              
         iii (3) = izz 
         jxx = iii (icut (1) ) 
         jyy = iii (icut (2) ) 
         izeig = offz (iz - 1) + (jxx - 1) * ny (iz) + jyy 
         z (izeig) = val (izz) 
      ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!------ fill X,Y arrays                                                 
!                                                                       
      DO i = 1, nx (iz) 
      x (offxy (iz - 1) + i) = xmin (iz) + (i - 1) * dxx 
      ENDDO 
      DO i = 1, ny (iz) 
      y (offxy (iz - 1) + i) = ymin (iz) + (i - 1) * dyy 
      ENDDO 
!                                                                       
!------ set other parameters                                            
!                                                                       
      lni (iz) = .true. 
      lh5 (iz) = .FALSE.
      lenc(iz) = max (nx (iz), ny (iz) ) 
      offxy (iz) = offxy (iz - 1) + lenc(iz) 
      offz (iz) = offz (iz - 1) + nx (iz) * ny (iz) 
      iz = iz + 1 
      CALL show_data (iz - 1) 
!                                                                       
      END SUBROUTINE read_4d                        
!*****7**************************************************************** 
      SUBROUTINE read_mess (ifil, ix, iy, str1, str2) 
!+                                                                      
!     Load MAN2 data files (MAN2 is a control program for the           
!     neutron diffractometers at FRM II. For more information           
!     contact the authors of KUPLOT.                                    
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!
      USE count_col_mod
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 20) 
!                                                                       
      CHARACTER ( * ) str1, str2 
      CHARACTER(LEN=PREC_STRING) :: line 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER ix, iy, ifil, nr, ianz, i, maxpp 
!                                                                       
!------ initial setup                                                   
!                                                                       
      achse (iwin, iframe, 1) = str1 
      achse (iwin, iframe, 2) = str2 
!                                                                       
!------ read data                                                       
!                                                                       
      nr = 1 
      maxpp = maxarray - offxy (iz - 1) 
      READ (ifil, 9999, end = 20) line 
      READ (ifil, 9999, end = 20) line 
   10 CONTINUE 
      READ (ifil, 9999, end = 20) line 
      CALL count_col (line, ianz) 
      IF (ianz.eq.0) goto 10 
      IF (ianz.gt.4) ianz = 4 
      BACKSPACE (ifil) 
   15 CONTINUE 
      READ (ifil, *, end = 20) (werte (i), i = 1, ianz) 
      IF (ix.gt.0) then 
         x (offxy (iz - 1) + nr) = werte (ix) 
      ELSE 
         x (offxy (iz - 1) + nr) = REAL(nr) 
      ENDIF 
      y (offxy (iz - 1) + nr) = werte (iy) 
      dx (offxy (iz - 1) + nr) = 0.0 
      dy (offxy (iz - 1) + nr) = 0.0 
      nr = nr + 1 
      IF (nr.gt.maxpp) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      GOTO 15 
   20 CONTINUE 
      lenc(iz) = nr - 1 
      offxy (iz) = offxy (iz - 1) + lenc(iz) 
      offz (iz) = offz (iz - 1) 
      iz = iz + 1 
!                                                                       
      CALL show_data (iz - 1) 
!                                                                       
 9999 FORMAT     (a) 
!                                                                       
      END SUBROUTINE read_mess                      
!*****7**************************************************************** 
      SUBROUTINE read_nipl (ifil, iwgb, istr, pgm) 
!+                                                                      
!     Read NIPL (KUPLOT own 2D format) and PGM (ASCII) files            
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER mrect
      PARAMETER (mrect = 50) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: line 
      CHARACTER(5) cmagic, kom 
      REAL rect (mrect, 4) 
      REAL dxx, dyy 
      INTEGER ifil, iwgb, istr, ir, i, j, zdummy 
      INTEGER maxpkt, maxzz 
      INTEGER nrect 
      LOGICAL pgm, inwg 
!                                                                       
!------ if there is a second file containing excluded regions,          
!------ read its contents                                               
!                                                                       
      IF (istr.eq.2.and..not.pgm) then 
         READ (iwgb, *, end = 50, err = 55) nrect 
         IF (nrect.gt.mrect) then 
            ier_num = - 9 
            ier_typ = ER_APPL 
            GOTO 999 
         ENDIF 
         DO i = 1, nrect 
         READ (iwgb, *, end = 50, err = 55) (rect (i, j), j = 1, 4) 
         ENDDO 
      ENDIF 
!                                                                       
      maxpkt = maxarray - offxy (iz - 1) 
      maxzz = maxarray - offz (iz - 1) 
!                                                                       
!------ read main data                                                  
!------ NIPL file format                                                
!                                                                       
      IF (.not.pgm) THEN
!
!        A header consisting of lines with a first '#' are ignored
         READ(ifil, '(a)', END=50, ERR=55) line
         DO WHILE(line(1:1)=='#')
            READ(ifil, '(a)', END=50, ERR=55) line
!           IF(line(1:18)=='# FOURIER : angles') READ(line(30:41),*) shear (1, 1)
!           IF(line(1:16)=='# FOURIER : aver'  ) THEN
!              READ(line(30:41),*) yskal_u (1, 1)
!              lyskal(1,1) = .TRUE.
!           ENDIF
         ENDDO
         READ (line, *, end = 50, err = 55) nx (iz), ny (iz) 
!        READ (ifil, *, end = 50, err = 55) nx (iz), ny (iz) 
!                                                                       
!------ - check array size                                              
!                                                                       
         IF (max (nx (iz), ny (iz) ) .gt.maxzz.or.max (nx (iz), ny (iz) &
         ) .gt.maxpkt) then                                             
            ier_num = - 6 
            ier_typ = ER_APPL 
            GOTO 999 
         ENDIF 
         READ (ifil, *, end = 50, err = 55) xmin (iz), xmax (iz),       &
         ymin (iz), ymax (iz)                                           
!                                                                       
!------ - check dimensions                                              
!                                                                       
         IF (xmin (iz) .eq.xmax (iz) .or.ymin (iz) .eq.ymax (iz) ) then 
            ier_num = - 29 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
         DO j = 1, ny (iz) 
         READ (ifil, *, end = 50, err = 55) (z (offz (iz - 1) + (i - 1) &
         * ny (iz) + j), i = 1, nx (iz) )                               
         ENDDO 
!                                                                       
!------ PGM file format                                                 
!                                                                       
      ELSE 
         READ (ifil, 1000, end = 50, err = 55) cmagic 
         IF (cmagic (1:2) .ne.'P2') then 
            ier_num = - 5 
            ier_typ = ER_APPL 
            GOTO 999 
         ENDIF 
!                                                                       
    5    READ (ifil, 1000, end = 50, err = 55) kom 
         IF (kom (1:1) .eq.'#') then 
            BACKSPACE (ifil) 
            READ (ifil, 1000, end = 50, err = 55) line 
            GOTO 5 
         ELSE 
            BACKSPACE (ifil) 
         ENDIF 
         READ (ifil, *, end = 50, err = 55) nx (iz), ny (iz) 
         READ (ifil, *, end = 50, err = 55) zdummy 
!                                                                       
!------ - check array size                                              
!                                                                       
         IF (max (nx (iz), ny (iz) ) .gt.maxzz.or.max (nx (iz), ny (iz) &
         ) .gt.maxpkt) then                                             
            ier_num = - 6 
            ier_typ = ER_APPL 
            GOTO 999 
         ENDIF 
         READ (ifil, *, end = 50, err = 55) ( (z (offz (iz - 1) +       &
         (i - 1) * ny (iz) + j), i = 1, nx (iz) ), j = ny (iz), 1,      &
         - 1)                                                           
         xmin (iz) = 0.0 
         xmax (iz) = REAL(nx (iz) - 1) 
         ymin (iz) = 0.0 
         ymax (iz) = REAL(ny (iz) - 1) 
      ENDIF 
!                                                                       
!------ set values for X and Y                                          
!                                                                       
      dxx = (xmax (iz) - xmin (iz) ) / REAL(nx (iz) - 1) 
      dyy = (ymax (iz) - ymin (iz) ) / REAL(ny (iz) - 1) 
      DO i = 1, nx (iz) 
      x (offxy (iz - 1) + i) = xmin (iz) + (i - 1) * dxx 
      ENDDO 
      DO i = 1, ny (iz) 
      y (offxy (iz - 1) + i) = ymin (iz) + (i - 1) * dyy 
      ENDDO 
!                                                                       
!------ check with excluded regions                                     
!                                                                       
      DO i = 1, nx (iz) 
      DO j = 1, ny (iz) 
      IF (istr.eq.2.and..not.pgm) then 
         inwg = .false. 
         DO ir = 1, nrect 
         inwg = inwg.or. (x (offxy (iz - 1) + i) .ge.rect (ir, 1)       &
         .and.x (offxy (iz - 1) + i) .le.rect (ir, 2) .and.y (offxy (iz &
         - 1) + j) .ge.rect (ir, 3) .and.y (offxy (iz - 1) + j)         &
         .le.rect (ir, 4) )                                             
         ENDDO 
         IF (inwg) z (offz (iz - 1) + (i - 1) * ny (iz) + j) = - 9999.0 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
!------ set remaining parameters                                        
!                                                                       
      lni (iz) = .true. 
      lh5 (iz) = .FALSE.
      lenc(iz) = max (nx (iz), ny (iz) ) 
      offxy (iz) = offxy (iz - 1) + lenc(iz) 
      offz (iz) = offz (iz - 1) + nx (iz) * ny (iz) 
      iz = iz + 1 
!                                                                       
      CALL show_data (iz - 1) 
      GOTO 999 
!                                                                       
   50 CONTINUE 
      ier_num = - 6 
      ier_typ = ER_IO 
      GOTO 999 
!                                                                       
   55 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
      GOTO 999 
!                                                                       
  999 CONTINUE 
 1000 FORMAT     (a) 
      END SUBROUTINE read_nipl                      
!*****7**************************************************************** 
      SUBROUTINE read_spec (ifil, ianz, cpara, lpara, werte, maxw,      &
      itype)                                                            
!+                                                                      
!     Load SPEC scan files                                              
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
USE precision_mod
USE str_comp_mod
!                                                                       
      IMPLICIT none 
!
      INTEGER, PARAMETER :: MAXWW = 1
!                                                                       
      INTEGER maxw, colm 
      INTEGER READ_TYPE_SC 
!     INTEGER READ_TYPE_ST 
      INTEGER READ_TYPE_SM 
!                                                                       
      PARAMETER (colm = 200) 
      PARAMETER (READ_TYPE_SC = 0) 
!     PARAMETER (READ_TYPE_ST = 1) 
      PARAMETER (READ_TYPE_SM = 2) 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER (LEN=PREC_STRING), DIMENSION(MAXWW) :: ccpara
      INTEGER             , DIMENSION(MAXWW) :: llpara
      REAL(KIND=PREC_DP)  , DIMENSION(MAXWW) :: wwerte
      CHARACTER(4096) mine, sinfo 
      CHARACTER(200) date, tit 
      CHARACTER(40) field (0:colm), input (colm) 
      CHARACTER(LEN=40) :: scan_type    ! Scan type on #S instruction
      CHARACTER(LEN=40) :: scan_mot     ! Motor type on #S instruction
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL col (0:colm), dummy 
      REAL mca_par (3) 
      REAL, DIMENSION(:,:), ALLOCATABLE :: scan_limits
      LOGICAL                           :: yes_limits
      INTEGER counts (8192) 
      INTEGER lpara (maxw), icell (5) 
      INTEGER ifil, ianz, iianz, il, ipt 
      INTEGER iscan, istart, iend, i, j 
      INTEGER nscan, nr, nf, maxpp 
      INTEGER nmca 
      INTEGER nscans 
      INTEGER itype 
      INTEGER  :: janz      ! temporary number of parameters
      INTEGER  :: npoints   ! number of data point on the #S instruction
      REAL     :: xstart    ! Start point          on the #S instruction
      REAL     :: xend      ! End   point          on the #S instruction
      REAL     :: ctime     ! Count time           on the #S instruction
      LOGICAL not_found, data_read, lsigma, lend, lall 
      LOGICAL lkev 
!                                                                       
      REAL chan2kev 
!                                                                       
      lend = .false. 
      yes_limits = .false.
      nr = 0
      lkev = .FALSE.
!                                                                       
!------ No further parameters -> list scans in file                     
!                                                                       
      nmca = 0 
      nscans = 0 
!                                                                       
      IF (ianz.eq.2) then 
         WRITE (output_io, 1000) fname (iz) (1:len_str (fname (iz) ) ) 
         ALLOCATE(scan_limits(2,1000))
    1    CONTINUE 
         READ (ifil, 9999, end = 2) mine 
         IF (mine (1:2) .eq.'#S') then 
            sinfo = mine 
            nscans = nscans + 1 
            ipt = 0 
   11       CONTINUE 
            READ (ifil, 9999, end = 22) mine 
            IF (mine (1:2) .ne.'#L') goto 11 
   12       CONTINUE 
            READ (ifil, 9999, end = 22) mine 
            IF (mine (1:2) .eq.'#C') then 
               GOTO 12 
            ELSEIF (mine (1:5) .eq.'#@MCA') then 
               CALL read_mca_point (ifil, mine, counts) 
               GOTO 12 
            ENDIF 
            READ (mine, *, end = 22, err = 22) dummy 
            ipt = ipt + 1 
            GOTO 12 
   22       CONTINUE 
            IF (ipt.gt.0) then 
               WRITE (output_io, 1100) ipt, sinfo (3:min (70, len_str ( &
               sinfo) ) )                                               
               CALL get_col (sinfo, field, nf, colm) 
               READ (field (1), * ) nscan 
               IF(nf==7) THEN
                  yes_limits = .true.
                  READ (field (2), * ) scan_type
                  READ (field (3), * ) scan_mot 
                  READ (field (4), * ) xstart 
                  READ (field (5), * ) xend 
                  READ (field (6), * ) npoints 
                  READ (field (7), * ) ctime
                  scan_limits(1,nscans) = xstart
                  scan_limits(2,nscans) = xend
               ENDIF 
               res_para (1) = nscans 
               res_para (  nscans + 2) = ipt 
               res_para (0) = 3*nscans + 2 
            ENDIF 
            BACKSPACE (ifil) 
         ELSEIF (mine (1:5) .eq.'#@MCA') then 
            CALL read_mca_point (ifil, mine, counts) 
            nmca = nmca + 1 
         ENDIF 
         GOTO 1 
    2    CONTINUE 
         IF (nmca.gt.0) then 
            WRITE (output_io, 1200) nmca 
         ENDIF 
         res_para (1) = nscans 
         res_para (2) = nmca 
         res_para (0) = nscans + 2 
         IF(yes_limits) THEN
         DO i=1, nscans
               res_para (1*nscans + 2 + i) = scan_limits(1,i)
               res_para (2*nscans + 2 + i) = scan_limits(2,i)
         ENDDO
         res_para (0) = 3*nscans + 2
         ENDIF 
         DEALLOCATE(scan_limits)
         RETURN 
!                                                                       
!------ Scan information                                                
!                                                                       
      ELSEIF (ianz.eq.4.and.str_comp (cpara (4) , 'info', 1, lpara (4) ,&
      4) ) then                                                         
         CALL read_spec_info (ifil, ianz, cpara, lpara, maxw) 
!       else                                                            
!         ier_num = -6                                                  
!         ier_typ = ER_COMM                                             
!         return                                                        
!       endif                                                           
!                                                                       
!------ Only scan number -> list scans in file                          
!                                                                       
      ELSEIF (ianz.eq.3.or. (ianz.ge.5.and.ianz.le.8) .or. (            &
      itype.eq.READ_TYPE_SM.and.ianz.eq.4) ) then                       
         data_read = ianz.gt.3 
         iianz = 1 
!                                                                       
!------ - Search for scan number first                                  
!                                                                       
         field (0) = 'ignored' 
         date = 'none' 
!                                                                       
         CALL get_scan_numbers (cpara (3), lpara (3), istart, iend,     &
         lall)                                                          
         IF (ier_num.ne.0) return 
!                                                                       
         iscan = istart 
!                                                                       
         CALL del_params (3, ianz, cpara, lpara, maxw) 
         DO i = 1, ianz 
         input (i) = cpara (i) (1:lpara (i) ) 
         ENDDO 
!                                                                       
 1111    CONTINUE 
!                                                                       
         IF (itype.eq.READ_TYPE_SM.and.istart.eq.0) then 
!                                                                       
!     ----Read a SPEC file that contains no scans, just MCA data        
!                                                                       
            CONTINUE 
         ELSE 
!                                                                       
!     ----Read a SPEC file that contains scans                          
!                                                                       
            not_found = .true. 
    3       CONTINUE 
            READ (ifil, 9999, end = 4) mine 
            IF (mine (1:2) .eq.'#S') then 
               CALL get_col (mine, field, nf, colm) 
               READ (field (1), * ) nscan 
               IF(nf==7) THEN
                  READ (field (2), * ) scan_type
                  READ (field (3), * ) scan_mot 
                  READ (field (4), * ) xstart 
                  READ (field (5), * ) xend 
                  READ (field (6), * ) npoints 
                  READ (field (7), * ) ctime
                  janz   = 1
                  ccpara(1) = input(1)
                  llpara(1) = LEN_TRIM(ccpara(1))
                  CALL ber_params (janz, ccpara, llpara, wwerte, MAXWW)
!                 linput = LEN_TRIM(input(1))
!                 CALL ber_params (janz, input, linput, werte, maxw) 
                  i = NINT(wwerte(1))
!                  READ (input(1),*) i
                  IF(nscan == iscan) THEN
                     IF(npoints+1 < i) THEN
                        ier_num = -41 
                        ier_typ = ER_APPL 
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF 
               IF (iscan.eq.nscan.or.lall) then 
                  tit = mine (1:MIN(200,LEN_TRIM(mine)))
   33             CONTINUE 
                  READ (ifil, 9999, end = 4) mine 
                  IF (mine (1:2) .eq.'#D') then 
                     date = mine (3:len_str (mine) ) 
                  ENDIF 
                  IF (mine (1:2) .eq.'#L') then 
                     CALL get_col (mine, field, nf, colm) 
                     not_found = .false. 
                     GOTO 4 
                  ENDIF 
                  GOTO 33 
               ENDIF 
            ENDIF 
            GOTO 3 
    4       CONTINUE 
            IF (not_found) then 
               ier_num = - 37 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
!                                                                       
!------ --- display info                                                
!                                                                       
            IF (date (1:4) .eq.'none') then 
               WRITE (output_io, 2000) iscan, fname (iz) (1:len_str (   &
               fname (iz) ) )                                           
            ELSE 
               WRITE (output_io, 2010) iscan, fname (iz) (1:len_str (   &
               fname (iz) ) ), date (1:len_str (date) )                 
            ENDIF 
!                                                                       
            IF (.not.data_read) then 
               WRITE (output_io, * ) 
               DO i = 1, nf 
               WRITE (output_io, 2020) i, field (i) (1:len_str (field ( &
               i) ) )                                                   
               ENDDO 
               WRITE (output_io, * ) 
               RETURN 
            ENDIF 
         ENDIF 
!                                                                       
         IF (itype.ne.READ_TYPE_SM) then 
            nr = 1
            DO j = 1, 5 
            icell (j) = 0 
            ENDDO 
!                                                                       
            DO i = 1, nf 
            il = len_str (field (i) ) 
            DO j = 1, ianz 
            IF (field (i) (1:il) .eq.input (j) (1:il) ) then 
               icell (j) = i 
               cpara (j) = '0' 
               lpara (j) = 1 
            ENDIF 
            ENDDO 
            ENDDO 
!                                                                       
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            not_found = .false. 
            DO j = 1, ianz 
            IF (icell (j) .eq.0) icell (j) = nint (werte (j) ) 
            IF (icell (j) .lt.0) icell (j) = nf + icell (j) + 1 
            not_found = not_found.or. (icell (j) .lt.0.or.icell (j)     &
            .gt.nf)                                                     
            ENDDO 
!                                                                       
            WRITE (output_io, 2100) 'x', field (icell (1) ) (1:len_str (&
            field (icell (1) ) ) ) , icell (1)                          
            WRITE (output_io, 2100) 'y', field (icell (2) ) (1:len_str (&
            field (icell (2) ) ) ) , icell (2)                          
            IF (icell (3) .ne.0) WRITE (output_io, 2110) field (icell ( &
            3) ) (1:len_str (field (icell (3) ) ) ), icell (3)          
            IF (icell (4) .ne.0) WRITE (output_io, 2120) 'dx', field (  &
            icell (4) ) (1:len_str (field (icell (4) ) ) ) , icell (4)  
            IF (icell (5) .ne.0) WRITE (output_io, 2120) 'dy', field (  &
            icell (5) ) (1:len_str (field (icell (5) ) ) ) , icell (5)  
            WRITE (output_io, * ) 
!                                                                       
            IF (not_found) then 
               ier_num = - 38 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
         ELSEIF (itype.eq.READ_TYPE_SM) then 
!                                                                       
!       read which point within the scan shall be extracted from MCA    
!                                                                       
            cpara (1) = input (1) 
            lpara (1) = len_str (cpara (1) ) 
            IF (str_comp (input (2) , 'kev', 3, len_str (input (2) ) ,  &
            3) ) then                                                   
               cpara (2) = input (3) 
               lpara (2) = len_str (input (3) ) 
               cpara (3) = input (4) 
               lpara (3) = len_str (input (4) ) 
               cpara (4) = input (5) 
               lpara (4) = len_str (input (5) ) 
               ianz = 4 
               lkev = .true. 
            ELSE 
               ianz = 1 
               lkev = .false. 
            ENDIF 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            icell (1) = nint (werte (1) ) 
            IF (lkev) then 
               mca_par (1) = werte (2) 
               mca_par (2) = werte (3) 
               mca_par (3) = werte (4) 
            ELSE 
               mca_par (1) = 0.0 
               mca_par (2) = 1.0 
               mca_par (3) = 0.0 
            ENDIF 
            nr = 0 
            READ (ifil, 9999, end = 20) mine 
            READ (ifil, 9999, end = 20) mine 
            DO 
            READ (ifil, 9999, end = 20) mine 
            IF (mine (1:5) .eq.'#@MCA') then 
               nr = nr + 1 
               IF (nr.eq.icell (1) ) then 
                  BACKSPACE (ifil) 
                  exit 
               ENDIF 
               ELSEIF(mine(1:2) == '#S') THEN
                  ier_num = -41 
                  ier_typ = ER_APPL 
                  RETURN
            ENDIF 
            ENDDO 
         ENDIF 
!                                                                       
!------ - Finally reading of the data                                   
!                                                                       
         col (0) = 0.0 
         lsigma = itype.eq.READ_TYPE_SC.and. (icell (5) .eq.0) 
         maxpp = maxarray - offxy (iz - 1) 
   10    CONTINUE 
         READ (ifil, 9999, end = 20) mine 
         IF (mine (1:2) .eq.'#C') then 
            GOTO 10 
         ELSEIF (mine (1:5) .eq.'#@MCA') then 
            CALL read_mca_point (ifil, mine, counts) 
            IF (itype.eq.READ_TYPE_SM) then 
               IF (nr.eq.icell (1) ) then 
                  DO i = 1, 8192 
                  x (offxy (iz - 1) + i) = chan2kev (i, mca_par) 
                  y (offxy (iz - 1) + i) = counts (i) 
                  dx (offxy (iz - 1) + i) = i 
                  dy (offxy (iz - 1) + i) = sqrt (real (counts (i) ) ) 
                  ENDDO 
                  nr = 8193 
                  GOTO 20 
               ENDIF 
            ENDIF 
            GOTO 10 
         ENDIF 
         IF (istart.ne.0) then 
            READ (mine, *, end = 20, err = 30) (col (i), i = 1, nf) 
         ENDIF 
         IF (itype.ne.READ_TYPE_SM) THEN 
            x (offxy (iz - 1) + nr) = col (icell (1) ) 
            y (offxy (iz - 1) + nr) = col (icell (2) ) 
            dx (offxy (iz - 1) + nr) = col (icell (4) ) 
            dy (offxy (iz - 1) + nr) = col (icell (5) ) 
            IF (lsigma.and.col (icell (2) ) .ge.0.0) then 
               dy (offxy (iz - 1) + nr) = sqrt (col (icell (2) ) ) 
            ENDIF 
            IF (icell (3) .ne.0.and.col (icell (3) ) .ne.0.0) then 
               y (offxy (iz - 1) + nr) = y (offxy (iz - 1) + nr)        &
               / col (icell (3) )                                       
               dy (offxy (iz - 1) + nr) = dy (offxy (iz - 1) + nr)      &
               / col (icell (3) )                                       
            ENDIF 
         ENDIF 
         nr = nr + 1 
         IF (nr.gt.maxpp) then 
            ier_num = - 6 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         GOTO 10 
   20    CONTINUE 
         lend = .true. 
   30    CONTINUE 
!                                                                       
!------ - Settings                                                      
!                                                                       
         lenc(iz) = nr - 1 
         offxy (iz) = offxy (iz - 1) + lenc(iz) 
         offz (iz) = offz (iz - 1) 
         iz = iz + 1 
!                                                                       
         IF (itype.ne.READ_TYPE_SM) then 
            achse (iwin, iframe, 1) = field (icell (1) ) (1:len_str (   &
            field (icell (1) ) ) )                                      
            IF (icell (3) .eq.0) then 
               achse (iwin, iframe, 2) = field (icell (2) ) (1:len_str (&
               field (icell (2) ) ) )                                   
            ELSE 
               achse (iwin, iframe, 2) = field (icell (2) ) (1:len_str (&
               field (icell (2) ) ) ) //'/'//field (icell (3) ) (1:     &
               len_str (field (icell (3) ) ) )                          
            ENDIF 
            titel (iwin, iframe, 1) (1:len_str (tit) ) = tit (1:len_str &
            (tit) )                                                     
            IF (date (1:4) .ne.'none') titel (iwin, iframe, 2) (1:      &
            len_str (date) ) = date (1:len_str (date) )                 
         ELSE 
!                                                                       
!     ----Settings for an MCA data point within the SPEC file           
!                                                                       
            DO i = 1, nf 
            res_para (i) = col (i) 
            ENDDO 
            res_para (0) = nf 
            IF (lkev) then 
               achse (iwin, iframe, 1) = 'keV' 
            ELSE 
               achse (iwin, iframe, 1) = 'Channel' 
            ENDIF 
            achse (iwin, iframe, 2) = 'Counts' 
            IF (istart.gt.0) then 
               titel (iwin, iframe, 1) (1:len_str (tit) ) = tit (1:     &
               len_str (tit) )                                          
            ELSE 
               WRITE (titel (iwin, iframe, 1), 6000) icell (1) 
            ENDIF 
            IF (date (1:4) .ne.'none') titel (iwin, iframe, 2) (1:      &
            len_str (date) ) = date (1:len_str (date) )                 
         ENDIF 
!                                                                       
         CALL show_data (iz - 1) 
!                                                                       
         IF (iscan.lt.iend.and..not.lend) then 
            BACKSPACE (ifil) 
!                                                                       
            iscan = iscan + 1 
            fname (iz) = fname (iz - 1) 
            GOTO 1111 
         ENDIF 
      ENDIF 
!                                                                       
 1000 FORMAT (' List of scans in file ',a,' :') 
 1100 FORMAT (i7,' pts --> #',a) 
 1200 FORMAT (' Number of MCA data points:',i7) 
 2000 FORMAT (1x,78('-'),/' Scan ',i4,' in file ',a,' :',/,1x,78('-')) 
 2010 FORMAT (1x,78('-'),/' Scan ',i4,' in file ',a,/,10x,              &
     &        ' measured ',a,/,1x,78('-'))                              
 2020 FORMAT ('   Data field ',i4,' : ',a) 
 2100 FORMAT ('   Axis ',a,' is       : ',a15,' (#',i4,')') 
 2110 FORMAT ('   Normalization   : ',a15,' (#',i4,')') 
 2120 FORMAT ('   Error column ',a,' : ',a15,' (#',i4,')') 
 6000 FORMAT ('MCA data point Nr. ',i7) 
 9999 FORMAT (a) 
!                                                                       
      END SUBROUTINE read_spec                      
!*****7**************************************************************** 
      SUBROUTINE read_spec_info (ifil, ianz, cpara, lpara, maxw) 
!+                                                                      
!     Read scan information from SPEC scan files                        
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw, maxr, maxc 
!                                                                       
      PARAMETER (maxr = 10) 
      PARAMETER (maxc = 15) 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER(4096) mine 
      CHARACTER(200) date, scan, energy 
      CHARACTER(40) field (0:maxc), field_t (0:maxr, maxc) 
      CHARACTER(40) out_t (5) 
      REAL field_v (0:maxr, maxc), out_v (5) 
      INTEGER lpara (maxw) 
      INTEGER ic, ir, nf, i, ii 
      INTEGER ifil, ianz, istart, iend, iscan, nscan 
      LOGICAL lall, lend 
!                                                                       
!                                                                       
      date = 'none' 
      energy = 'none' 
!                                                                       
      DO ir = 0, maxr 
      DO ic = 1, maxc 
      field_t (ir, ic) = 'none' 
      field_v (ir, ic) = 0.0 
      ENDDO 
      ENDDO 
!                                                                       
!------ Get header information                                          
!                                                                       
    1 CONTINUE 
      READ (ifil, 9999, end = 2) mine 
      IF (mine (1:2) .eq.'#S') then 
         BACKSPACE (ifil) 
         GOTO 2 
      ELSEIF (mine (1:2) .eq.'#D') then 
         date = mine (3:len_str (mine) ) 
      ELSEIF (mine (1:6) .eq.'#C E =') then 
         energy = mine (3:len_str (mine) ) 
      ELSEIF (mine (1:2) .eq.'#O') then 
         READ (mine (3:4), * ) ir 
         CALL get_col (mine, field, nf, maxc) 
         DO ic = 1, nf 
         field_t (ir, ic) = field (ic) 
         ENDDO 
      ENDIF 
      GOTO 1 
    2 CONTINUE 
!                                                                       
!------ Search for scan number first                                    
!                                                                       
      CALL get_scan_numbers (cpara (3), lpara (3), istart, iend, lall) 
      IF (ier_num.ne.0) return 
!                                                                       
      iscan = istart 
!                                                                       
 1111 CONTINUE 
      lend = .true. 
!                                                                       
    3 CONTINUE 
      READ (ifil, 9999, end = 20) mine 
      IF (mine (1:2) .eq.'#S') then 
         scan = mine (3:len_str (mine) ) 
         CALL get_col (mine, field, nf, maxc) 
         READ (field (1), * ) nscan 
         lend = .false. 
         IF (iscan.eq.nscan.or.lall) goto 4 
      ENDIF 
      GOTO 3 
    4 CONTINUE 
!                                                                       
!------ display info                                                    
!                                                                       
      WRITE (output_io, 2000) iscan, fname (iz) (1:len_str (fname (iz) )&
      )                                                                 
!                                                                       
!------ read information                                                
!                                                                       
   10 CONTINUE 
      READ (ifil, 9999, end = 20) mine 
      IF (mine (1:1) .ne.'#') then 
         GOTO 30 
      ELSEIF (mine (1:2) .eq.'#D') then 
         date = mine (3:len_str (mine) ) 
      ELSEIF (mine (1:2) .eq.'#P') then 
         READ (mine (3:4), * ) ir 
         CALL get_col (mine, field, nf, maxc) 
         DO ic = 1, nf 
         READ (field (ic), * ) field_v (ir, ic) 
         ENDDO 
      ENDIF 
      GOTO 10 
   20 CONTINUE 
      lend = .true. 
   30 CONTINUE 
!                                                                       
      IF (lend) return 
!                                                                       
!------ display result                                                  
!                                                                       
      WRITE (output_io, 3000) scan (1:len_str (scan) ) 
      IF (date (1:4) .ne.'none') WRITE (output_io, 3010) date (1:       &
      len_str (date) )                                                  
      IF (energy (1:4) .ne.'none') WRITE (output_io, 3020) energy (1:   &
      len_str (energy) )                                                
!                                                                       
      ii = 0 
!                                                                       
      WRITE (output_io, 4000) 
      DO ir = 0, maxr 
      DO ic = 1, maxc 
      IF (field_t (ir, ic) .ne.'none') then 
         ii = ii + 1 
         out_t (ii) = field_t (ir, ic) 
         out_v (ii) = field_v (ir, ic) 
      ENDIF 
      IF (ii.eq.5) then 
         WRITE (output_io, 4010) (out_t (i), i = 1, 5) 
         WRITE (output_io, 4020) (out_v (i), i = 1, 5) 
         ii = 0 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      IF (ii.ne.0) then 
         WRITE (output_io, 4010) (out_t (i), i = 1, ii) 
         WRITE (output_io, 4020) (out_v (i), i = 1, ii) 
      ENDIF 
!                                                                       
!------ goto next scan                                                  
!                                                                       
      IF (iscan.lt.iend.and..not.lend) then 
         BACKSPACE (ifil) 
         iscan = iscan + 1 
         GOTO 1111 
      ENDIF 
!                                                                       
 2000 FORMAT (1x,78('-'),/' Scan ',i4,' in file ',a,' :',/,1x,78('-')) 
 3000 FORMAT (3x,'Scan       : ',a) 
 3010 FORMAT (3x,'Date       : ',a) 
 3020 FORMAT (3x,'Energy     : ',a) 
 4000 FORMAT (3x,'Parameters : ',/) 
 4010 FORMAT (3x,5(a13,3x)) 
 4020 FORMAT (3x,5(g13.6,3x)) 
 9999 FORMAT (a) 
!                                                                       
      END SUBROUTINE read_spec_info                 
!*****7**************************************************************** 
      SUBROUTINE get_scan_numbers (str, length, istart, iend, lall) 
!+                                                                      
!     Analyse scan string (number or number>number)                     
!-                                                                      
      USE ber_params_mod
USE lib_length
USE precision_mod
USE str_comp_mod
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) str 
      INTEGER length, istart, iend 
      LOGICAL lall 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ipos
!                                                                       
!                                                                       
      IF (str_comp (str, 'all', 1, length, 3) ) then 
         istart = 1 
         iend = 9999 
         lall = .true. 
      ELSE 
         ipos = index (str (1:length) , '>') 
         IF (ipos.eq.0) then 
            cpara (1) = str 
            lpara (1) = length 
            cpara (2) = str 
            lpara (2) = length 
         ELSE 
            cpara (1) = str (1:ipos - 1) 
            lpara (1) = len_str (cpara (1) ) 
            cpara (2) = str (ipos + 1:length) 
            lpara (2) = len_str (cpara (2) ) 
         ENDIF 
!                                                                       
         CALL ber_params (2, cpara, lpara, werte, maxw) 
         istart = nint (werte (1) ) 
         iend = nint (werte (2) ) 
         lall = .false. 
      ENDIF 
!                                                                       
      END SUBROUTINE get_scan_numbers               
!*****7**************************************************************** 
      SUBROUTINE read_mca_point (ifil, mine, counts) 
!+                                                                      
!     Reads a single MCA data point                                     
!-                                                                      
USE lib_length
      IMPLICIT none 
!                                                                       
      CHARACTER(4096) mine 
      INTEGER ifil 
      INTEGER counts (8192) 
!                                                                       
      INTEGER nchan 
      INTEGER chana 
      INTEGER chane 
      INTEGER nlines 
      INTEGER nrest 
      INTEGER ibsl 
      INTEGER :: length 
      INTEGER nl, i, j 
      INTEGER  :: iostatus
!                                                                       
!                                                                       
      DO i = 1, 8192 
      counts (i) = 0 
      ENDDO 
!                                                                       
      READ (ifil, 9999, end = 20) mine 
      READ (mine (8:len_str (mine) ), * ) nchan, chana, chane 
!                                                                       
      i = chana 
      READ (ifil, 9999, end = 20) mine 
!                                                                       
      DO while (mine (1:1) .eq.'#') 
      READ (ifil, 9999, end = 20) mine 
      ENDDO 
      ibsl = index (mine, '\') 
      READ (mine (3:ibsl - 1), * , IOSTAT=iostatus) (counts (j), j = i, i + 15) 
      IF(IS_IOSTAT_END(iostatus) ) GOTO 20
      IF(IS_IOSTAT_EOR(iostatus) ) GOTO 20
      nlines = (chane-chana + 1) / 16 - 1 
      nrest = mod (chane-chana, 16) + 1 
!                                                                       
      DO nl = 1, nlines 
      READ (ifil, 9999, end = 20) mine 
         length = LEN(TRIM(mine))
      i = i + 16 
      ibsl = index (mine, '\') 
      IF(ibsl == 0) ibsl = length+1
      READ (mine (1:ibsl - 1), * , IOSTAT=iostatus) (counts (j), j = i, i + 15) 
      IF(IS_IOSTAT_END(iostatus) ) GOTO 20
      IF(IS_IOSTAT_EOR(iostatus) ) GOTO 20
      ENDDO 
      IF(0 < nrest .AND. nrest < 16) THEN
         i = i + 16 
         READ (ifil, 9999, end = 20) mine 
         length = LEN(TRIM(mine))
         ibsl = index (mine, '\') 
         IF(ibsl == 0) ibsl = length+1
         READ (mine (1:ibsl - 1), * ) (counts (j), j = i, i + nrest - 1) 
      ENDIF
!                                                                       
   20 CONTINUE 
      RETURN 
 9999 FORMAT    (a) 
!                                                                       
      END SUBROUTINE read_mca_point                 
!*****7**************************************************************** 
      SUBROUTINE read_gsas (ifil, ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     Load GSAS files                                                   
!-                                                                      
      USE build_name_mod
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
USE precision_mod
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      INTEGER nbank 
!                                                                       
      PARAMETER (nbank = 100) 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER(LEN=PREC_STRING) :: line 
      CHARACTER(LEN=PREC_STRING) :: iname, gsastit, filename 
      CHARACTER(1) units, cunits 
      REAL(8) difc (nbank) 
      REAL(8) difa (nbank) 
      REAL(8) tzero (nbank) 
      REAL twotheta (nbank) 
      REAL vrange (nbank, 2), vp (nbank, 12) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL xval (maxarray) 
      REAL yval (maxarray) 
      REAL dyval (maxarray) 
      INTEGER lpara (maxw) 
      INTEGER bid (nbank) 
      INTEGER ifil, ianz 
      INTEGER iscan, istart, iend, i 
      INTEGER ndat, maxpp 
      INTEGER vtype (nbank) 
      LOGICAL lend, lall, liparm, lunits, lnorm 
!                                                                       
      INTEGER gsas_no_banks 
!                                                                       
      filename = fname (iz) 
!                                                                       
      lend = .false. 
      liparm = .false. 
      lunits = .false. 
!                                                                       
      DO i = 1, nbank 
      bid (i) = i 
      twotheta (i) = - 9999. 
      ENDDO 
!                                                                       
      IF (ianz.eq.2) then 
!                                                                       
!------ - Just print bank headers                                       
!                                                                       
   55    CONTINUE 
         READ (ifil, '(a)', end = 60) line 
         IF (line (1:4) .eq.'BANK') then 
            WRITE (output_io, 5000) line (1:len_str (line) ) 
         ENDIF 
         GOTO 55 
   60    CONTINUE 
      ELSEIF (ianz.ge.3) then 
!                                                                       
!------ - First get desired bank numbers                                
!                                                                       
         CALL get_scan_numbers (cpara (3), lpara (3), istart, iend,     &
         lall)                                                          
         IF (ier_num.ne.0) return 
         IF (lall) iend = gsas_no_banks (ifil, bid, nbank) 
!                                                                       
!------ - Now check units parameter                                     
!                                                                       
         IF (ianz.gt.3) then 
            units = cpara (4) (1:1) 
            CALL do_cap (units) 
            lunits = .true. 
         ENDIF 
!                                                                       
!------ - Now check for instrument parameter file and read              
!                                                                       
         IF (ianz.gt.4) then 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 5) 
            IF (ier_num.ne.0) return 
            iname = cpara (5) 
            WRITE (output_io, 1000) iname (1:len_str (iname) ) 
            CALL read_iparm (iname, nbank, difc, difa, tzero, twotheta, &
            vtype, vrange, vp)                                          
            IF (ier_num.ne.0) return 
            liparm = .true. 
!                                                                       
         ELSEIF (ianz.eq.4) then 
            CALL get_iparm_name (ifil, iname, liparm) 
            IF (liparm) then 
               WRITE (output_io, 1000) iname (1:len_str (iname) ) 
               CALL read_iparm (iname, nbank, difc, difa, tzero,        &
               twotheta, vtype, vrange, vp)                             
               IF (ier_num.ne.0) return 
            ENDIF 
         ENDIF 
!                                                                       
!------ - Check form norm flag                                          
!                                                                       
         lnorm = .false. 
         IF (ianz.gt.5) then 
            CALL do_cap (cpara (6) ) 
            lnorm = (cpara (6) (1:1) .eq.'N') 
         ENDIF 
!                                                                       
!------ - Now Read the data                                             
!                                                                       
         iscan = istart 
!                                                                       
 1111    CONTINUE 
         CALL read_gsas_data (ifil, iscan, xval, yval, dyval, ndat, bid,&
         nbank, cunits, gsastit)                                        
         IF (ier_num.ne.0) return 
!                                                                       
         CALL gsas_norm (xval, yval, dyval, ndat, iscan, bid, vtype,    &
         vrange, vp, lnorm, nbank)                                      
         IF (ier_num.ne.0) return 
!                                                                       
         CALL convert_units (xval, yval, dyval, ndat, units, cunits,    &
         lunits, difc (bid (iscan) ), difa (bid (iscan) ), tzero (bid ( &
         iscan) ), twotheta (bid (iscan) ) )                            
         IF (ier_num.ne.0) return 
!                                                                       
         maxpp = maxarray - offxy (iz - 1) 
         IF (ndat.gt.maxpp) then 
            ier_num = - 6 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
         DO i = 1, ndat 
         x (offxy (iz - 1) + i) = xval (i) 
         y (offxy (iz - 1) + i) = yval (i) 
         dx (offxy (iz - 1) + i) = 0.0 
         dy (offxy (iz - 1) + i) = dyval (i) 
         ENDDO 
!                                                                       
!------ - Settings                                                      
!                                                                       
         lenc(iz) = ndat 
         offxy (iz) = offxy (iz - 1) + lenc(iz) 
         offz (iz) = offz (iz - 1) 
         fform (iz) = 'XY' 
         iz = iz + 1 
!                                                                       
      IF (units.eq.'T') achse (iwin, iframe, 1)  = 'Time of flight (\\gm&
     &sec)'                                                             
         IF (units.eq.'Q') achse (iwin, iframe, 1) = 'Q (\\A\\u-1\\d)' 
         IF (units.eq.'D') achse (iwin, iframe, 1) = 'd spacing (\\A)' 
         IF (units.eq.'L') achse (iwin, iframe, 1) = '\\gl (\\A)' 
         IF (units.eq.'X') achse (iwin, iframe, 1) = '2\\gH (degrees)' 
         achse (iwin, iframe, 2) = 'Intensity' 
         titel (iwin, iframe, 1) = gsastit (1:len_str (gsastit) ) 
         IF (bid (iscan) .lt.10) then 
            WRITE (fname (iz - 1), 2000) filename (1:len_str (filename) &
            ), bid (iscan)                                              
         ELSE 
            WRITE (fname (iz - 1), 2100) filename (1:len_str (filename) &
            ), bid (iscan)                                              
         ENDIF 
!                                                                       
         CALL show_data (iz - 1) 
!                                                                       
         IF (iscan.lt.iend.and..not.lend) then 
            BACKSPACE (ifil) 
            iscan = iscan + 1 
            fname (iz) = filename (1:MIN(200,LEN_TRIM(filename)))
            GOTO 1111 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (' Reading GSAS instrument parameter file ',a, ' ..') 
 2000 FORMAT     (a,'.Bank_',i1) 
 2100 FORMAT     (a,'.Bank_',i2) 
 5000 FORMAT     (1x,a) 
      END SUBROUTINE read_gsas                      
!*****7**************************************************************** 
      SUBROUTINE gsas_norm (xval, yval, dyval, ndat, ib, bid, vtype,    &
      vrange, vp, lnor, nbank)                                          
!+                                                                      
!     Normalizes by incident spectrum                                   
!-                                                                      
      USE debug_mod 
      USE errlist_mod 
      USE wink_mod
!
      USE prompt_mod 
      USE kuplot_config 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER nbank 
!                                                                       
      REAL xval (maxarray) 
      REAL yval (maxarray) 
      REAL dyval (maxarray) 
      REAL vrange (nbank, 2), vp (nbank, 12) 
      REAL cof (12) 
      REAL i0, tof 
      INTEGER vtype (nbank) 
      INTEGER bid (nbank) 
      INTEGER ndat, i, j, ib, ii 
      LOGICAL lnor 
!                                                                       
      IF (.not.lnor) return 
!                                                                       
      IF (vtype (bid (ib) ) .lt.1.or.vtype (bid (ib) ) .gt.4) then 
         ier_typ = ER_APPL 
         ier_num = - 52 
         RETURN 
      ENDIF 
!                                                                       
      WRITE (output_io, 1000) 
      IF (dbg) open (91, file = "i0.dat", status = 'unknown') 
!                                                                       
!------ Cut range to vrange                                             
!                                                                       
      ii = 1 
      DO i = 1, ndat 
      tof = xval (i) / 1000. 
      IF (tof.gt.vrange (bid (ib), 1) .and.tof.lt.vrange (bid (ib),     &
      2) ) then                                                         
         xval (ii) = xval (i) 
         yval (ii) = yval (i) 
         dyval (ii) = dyval (i) 
         ii = ii + 1 
      ENDIF 
      ENDDO 
!                                                                       
      ndat = ii - 1 
!                                                                       
      DO i = 1, ndat 
      tof = xval (i) / 1000. 
      IF (vtype (bid (ib) ) .eq.1) then 
         i0 = vp (bid (ib), 1) 
         i0 = i0 + vp (bid (ib), 2) * exp ( - 1.0 * vp (bid (ib),       &
         3) * tof)                                                      
         i0 = i0 + vp (bid (ib), 4) * exp ( - 1.0 * vp (bid (ib),       &
         5) * tof**2)                                                   
         i0 = i0 + vp (bid (ib), 6) * exp ( - 1.0 * vp (bid (ib),       &
         7) * tof**3)                                                   
         i0 = i0 + vp (bid (ib), 8) * exp ( - 1.0 * vp (bid (ib),       &
         9) * tof**4)                                                   
         i0 = i0 + vp (bid (ib), 10) * exp ( - 1.0 * vp (bid (ib),      &
         11) * tof**5)                                                  
!                                                                       
      ELSEIF (vtype (bid (ib) ) .eq.2) then 
         i0 = vp (bid (ib), 1) 
         i0 = i0 + (vp (bid (ib), 2) / tof**5) * exp ( - 1.0 * vp (bid (&
         ib), 3) / tof**2)                                              
         i0 = i0 + vp (bid (ib), 4) * exp ( - 1.0 * vp (bid (ib),       &
         5) * tof**2)                                                   
         i0 = i0 + vp (bid (ib), 6) * exp ( - 1.0 * vp (bid (ib),       &
         7) * tof**3)                                                   
         i0 = i0 + vp (bid (ib), 8) * exp ( - 1.0 * vp (bid (ib),       &
         9) * tof**4)                                                   
         i0 = i0 + vp (bid (ib), 10) * exp ( - 1.0 * vp (bid (ib),      &
         11) * tof**5)                                                  
!                                                                       
      ELSEIF (vtype (bid (ib) ) .eq.3) then 
         CALL calc_cheb (tof, cof) 
         i0 = vp (bid (ib), 1) 
         DO j = 2, 12 
         i0 = i0 + vp (bid (ib), j) * cof (j) 
         ENDDO 
!                                                                       
      ELSEIF (vtype (bid (ib) ) .eq.4) then 
         CALL calc_cheb (tof, cof) 
         i0 = vp (bid (ib), 1) 
         i0 = i0 + (vp (bid (ib), 2) / tof**5) * exp ( - 1.0 * vp (bid (&
         ib), 3) / tof**2)                                              
         DO j = 2, 10 
         i0 = i0 + vp (bid (ib), j + 2) * cof (j) 
         ENDDO 
      ENDIF 
      yval (i) = yval (i) / i0 
      dyval (i) = dyval (i) / i0 
!                                                                       
      IF (dbg) WRITE (91, * ) xval (i), i0 
      ENDDO 
!                                                                       
      IF (dbg) close (91) 
!                                                                       
 1000 FORMAT     (' Normalizing by incident spectrum ..') 
      END SUBROUTINE gsas_norm                      
!*****7**************************************************************** 
      SUBROUTINE calc_cheb (tof, cof) 
!+                                                                      
!     Calculated Chebyschev polynominal coefficients                    
!-                                                                      
      IMPLICIT none 
!                                                                       
      REAL tof, cof (12) 
      REAL t 
      INTEGER i 
!                                                                       
      t = 1.0 / tof 
      t = 2.0 * t - 1.0 
!                                                                       
      cof (1) = 1.0 
      cof (2) = t 
!                                                                       
      DO i = 3, 12 
      cof (i) = 2.0 * t * cof (i - 1) - cof (i - 2) 
      ENDDO 
!                                                                       
      END SUBROUTINE calc_cheb                      
!*****7**************************************************************** 
      SUBROUTINE convert_units (xval, yval, dyval, ndat, units, cunits, &
      lunits, difc, difa, tzero, tth)                                   
!+                                                                      
!     Convert units from TOF to D,L,Q                                   
!     For X we convert T to degrees                                     
!-                                                                      
      USE debug_mod 
      USE errlist_mod 
USE lib_length
      USE wink_mod
      USE prompt_mod 
      USE kuplot_config 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(1) units, cunits 
      REAL(8) difc, difa, tzero 
      REAL(8) secondterm 
      REAL xval (maxarray) 
      REAL yval (maxarray) 
      REAL dyval (maxarray) 
      REAL tth, t, d, tsint 
      INTEGER ndat, j 
      LOGICAL lunits 
!                                                                       
!                                                                       
!------ Keep units                                                      
!                                                                       
      IF (.not.lunits) then 
         units = cunits 
         RETURN 
      ENDIF 
!                                                                       
!------ Selected same units                                             
!                                                                       
      IF (units.eq.cunits) return 
!                                                                       
!------ --------------------------------------------------------------- 
!------ Go from TOF to DEGEES (it is XRAY data !)                       
!------ --------------------------------------------------------------- 
!                                                                       
      IF (cunits.eq.'T'.and.units.eq.'X') then 
         WRITE (output_io, 1100) units (1:len_str (units) ) 
         DO j = 1, ndat 
         xval (j) = 0.01 * xval (j) 
         ENDDO 
         RETURN 
      ENDIF 
!                                                                       
!------ --------------------------------------------------------------- 
!------ Go from TOF to D,Q,L                                            
!------ --------------------------------------------------------------- 
!                                                                       
      IF (cunits.eq.'T') then 
         WRITE (output_io, 1100) units (1:len_str (units) ) 
         IF (tth.eq. - 9999.) then 
            WRITE (output_io, 1200) 
            units = 'T' 
            RETURN 
         ENDIF 
!                                                                       
!------ - D spacing                                                     
!------ - Use d = [-DIFC +/- sqrt(DIFC**2 + 4*DIFA*(TOF-TZERO))]/2*DIFA 
!------ -     In case DIFA is 0.0 use d=(T-T0)/DIFC                     
!                                                                       
         IF (units.eq.'D') then 
            IF (difa.ne.0.0) then 
               DO j = 1, ndat 
               t = REAL(REAL(xval (j),KIND=KIND(0.0D0)) - tzero) 
               secondterm = dsqrt (4.0 * difa * t + difc * difc) 
               xval (j) = REAL(( - difc + secondterm) / 2.D0 / difa)
               ENDDO 
            ELSE 
               DO j = 1, ndat 
               xval (j) = REAL((REAL(xval (j),KIND=KIND(0.0D0)) - tzero) / difc )
               ENDDO 
            ENDIF 
!                                                                       
!------ - Wavelength                                                    
!------ - Use l = 2d sin(theta)                                         
!                                                                       
         ELSEIF (units.eq.'L') then 
            tsint = 2. * sin (abs (tth * REAL(pi)) / 360.) 
            IF (difa.ne.0.0) then 
               DO j = 1, ndat 
               t = REAL(REAL(xval (j),KIND=KIND(0.0D0)) - tzero) 
               secondterm = dsqrt (4.0D0 * difa * REAL(t,KIND=KIND(0.0D0)) + difc * difc) 
               d = REAL(( - difc + secondterm) / 2.D0 / difa) 
               xval (j) = d * tsint 
               ENDDO 
            ELSE 
               DO j = 1, ndat 
               d = REAL((REAL(xval (j),KIND=KIND(0.0D0)) - tzero) / difc) 
               xval (j) = d * tsint 
               ENDDO 
            ENDIF 
!                                                                       
!------ - Q                                                             
!------ - Use Q = 2 pi / d (and reverse arrays)                         
!                                                                       
         ELSEIF (units.eq.'Q') then 
            IF (difa.ne.0.0) then 
               DO j = 1, ndat 
               t = REAL(REAL(xval (j),KIND=KIND(0.0D0)) - tzero) 
               secondterm = dsqrt (4.0 * difa * REAL(t,KIND=KIND(0.0D0)) + difc * difc) 
               d = REAL(( - difc + secondterm) / 2.D0 / difa )
               xval (j) = 2 * REAL(pi) / d 
               ENDDO 
            ELSE 
               DO j = 1, ndat 
               d = REAL((REAL(xval (j),KIND=KIND(0.0D0)) - tzero) / difc) 
               xval (j) = 2 * REAL(pi) / d 
               ENDDO 
            ENDIF 
!                                                                       
            CALL reverse_array (xval, ndat) 
            CALL reverse_array (yval, ndat) 
            CALL reverse_array (dyval, ndat) 
!                                                                       
!------ - No valid unit                                                 
!                                                                       
         ELSE 
            ier_num = - 51 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ --------------------------------------------------------------- 
!------ Go from D to Q                                                  
!------ --------------------------------------------------------------- 
!                                                                       
      ELSEIF (cunits.eq.'D') then 
         WRITE (output_io, 1100) units (1:len_str (units) ) 
!                                                                       
!------ - D to Q using Q=2 pi / D                                       
!                                                                       
         IF (units.eq.'Q') then 
            DO j = 1, ndat 
            xval (j) = 2 * REAL(pi) / xval (j) 
            ENDDO 
!                                                                       
            CALL reverse_array (xval, ndat) 
            CALL reverse_array (yval, ndat) 
            CALL reverse_array (dyval, ndat) 
!                                                                       
!------ - No valid unit                                                 
!                                                                       
         ELSE 
            ier_num = - 51 
            ier_typ = ER_APPL 
         ENDIF 
!------ --------------------------------------------------------------- 
!------ Go from Q to D                                                  
!------ --------------------------------------------------------------- 
!                                                                       
      ELSEIF (cunits.eq.'Q') then 
         WRITE (output_io, 1100) units (1:len_str (units) ) 
!                                                                       
!------ - D to Q using Q=2 pi / D                                       
!                                                                       
         IF (units.eq.'D') then 
            DO j = 1, ndat 
            xval (j) = 2 * REAL(pi) / xval (j) 
            ENDDO 
!                                                                       
            CALL reverse_array (xval, ndat) 
            CALL reverse_array (yval, ndat) 
            CALL reverse_array (dyval, ndat) 
!                                                                       
!------ - No valid unit                                                 
!                                                                       
         ELSE 
            ier_num = - 51 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
!                                                                       
 1100 FORMAT     (' Converting to units ',a, ' ..') 
 1200 FORMAT     (' Unit conversion requires instrument parameter',     &
     &                   ' file ..')                                    
      END SUBROUTINE convert_units                  
!*****7**************************************************************** 
      SUBROUTINE reverse_array (dat, ndat) 
!+                                                                      
!     Reverse given array                                               
!-                                                                      
      USE kuplot_config 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL dat (maxarray) 
      REAL dummy (maxarray) 
      INTEGER i, ndat 
!                                                                       
      DO i = 1, ndat 
      dummy (i) = dat (i) 
      ENDDO 
!                                                                       
      DO i = 1, ndat 
      dat (i) = dummy (ndat - i + 1) 
      ENDDO 
!                                                                       
      END SUBROUTINE reverse_array                  
!*****7**************************************************************** 
      SUBROUTINE get_iparm_name (ifil, iname, liparm) 
!+                                                                      
!     Read GSAS data                                                    
!-                                                                      
      USE debug_mod 
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
USE lib_length
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) iname 
      INTEGER ifil 
      LOGICAL liparm 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: line, dummy 
!                                                                       
      liparm = .false. 
      iname = '' 
!                                                                       
      REWIND (ifil) 
!                                                                       
!----- Search for Instrument parameter file line                        
!                                                                       
   10 CONTINUE 
      READ (ifil, '(a)', end = 15) line 
      dummy = line 
      CALL do_cap (dummy) 
      IF (dummy (1:26) .eq.'INSTRUMENT PARAMETER FILE:') then 
         iname = line (27:len_str (dummy) ) 
         liparm = .true. 
      ENDIF 
      IF (.not.liparm) goto 10 
   15 CONTINUE 
!                                                                       
!----- Check the file exists                                            
!                                                                       
      INQUIRE (file = iname, exist = liparm) 
      IF (.not.liparm) WRITE (output_io, 1000) iname (1:len_str (iname) &
      )                                                                 
!                                                                       
 1000 FORMAT (' Instrument parameter file ',a,' not found ..') 
      END SUBROUTINE get_iparm_name                 
!*****7**************************************************************** 
      SUBROUTINE read_gsas_data (ifil, ibank, xval, yval, dyval, ndat,  &
      bid, nbank, cunits, gsastit)                                      
!+                                                                      
!     Read GSAS data                                                    
!-                                                                      
      USE debug_mod 
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER nbank 
!                                                                       
      CHARACTER ( * ) gsastit 
      CHARACTER(1) cunits 
      REAL xval (maxarray) 
      REAL yval (maxarray) 
      REAL dyval (maxarray) 
      INTEGER bid (nbank) 
      INTEGER ndat, ibank, ifil 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: line, bank_line 
      REAL clckwdt, onetick, scal 
      REAL c_min, c_d, monitor 
      REAL gtmin, gtmax, gtlog, xdelta 
      INTEGER tmap (3, maxarray), tmap3i, ntmax 
      INTEGER i, j, ll, ib, nc, nr, iq, lt 
      INTEGER idummy, itmap, icons, ilog 
      INTEGER itmap1i, itmap2i, itmap3i 
      INTEGER itmap1iplus1 
      INTEGER i_tmap 
      INTEGER maxmapu 
      INTEGER it_no, it_nval, it_nrec, imax 
      LOGICAL ltmap, lbank, lmon 
!                                                                       
!                                                                       
      ltmap = .false. 
      lbank = .false. 
      lmon = .false. 
!                                                                       
      REWIND (ifil) 
!                                                                       
!----- Read title line                                                  
!                                                                       
      READ (ifil, '(a)', end = 15, err = 998) gsastit 
!                                                                       
!----- Search for Monitor count                                         
!                                                                       
   10 CONTINUE 
      READ (ifil, '(a)', end = 15, err = 998) line 
      ll = len_str (line) 
      IF (line (1:8) .eq.'Monitor:'.or.line (1:8) .eq.'MONITOR:') then 
         READ (line (9:ll), *, err = 998) monitor 
         lmon = .true. 
      ENDIF 
      IF (.not.lmon) goto 10 
                                                                        
      IF (lmon) then 
         res_para (0) = 1 
         res_para (1) = monitor 
         WRITE (output_io, 5000) monitor 
      ENDIF 
!                                                                       
!----- Now try to find the bank we are looking for                      
!                                                                       
   15 CONTINUE 
      REWIND (ifil) 
      scal = 1.0 
   20 CONTINUE 
      READ (ifil, '(a)', end = 510, err = 998) line 
      ll = len_str (line) 
      IF (line (1:15) .eq.'#scale_factor =') then 
         READ (line (16:ll), *, err = 998) scal 
      ENDIF 
      IF (line (1:4) .eq.'BANK') then 
         itmap = index (line (5:ll) , 'TIME_MAP') 
         icons = index (line (5:ll) , 'CON') 
         ilog = index (line (5:ll) , 'SLOG') 
!                                                                       
         IF (index (line (5:ll) , 'TIME_MAP') .ne.0) cunits = 'T' 
         IF (index (line (5:ll) , 'SLOG') .ne.0) cunits = 'T' 
         IF (index (line (5:ll) , 'CONS') .ne.0) cunits = 'T' 
         IF (index (line (5:ll) , 'CONQ') .ne.0) cunits = 'Q' 
         IF (index (line (5:ll) , 'COND') .ne.0) cunits = 'D' 
!                                                                       
!----- --- Type TIME_MAP                                                
                                                                        
         IF (itmap.gt.0) then 
            READ (line (5:itmap + 3), *, err = 998) ib, nc, nr 
            lbank = (ib.eq.bid (ibank) ) 
            IF (lbank) then 
               bank_line = line 
               READ (line (itmap + 12:ll), *, err = 998) i_tmap 
               IF (dbg) then 
                  WRITE ( *, 1000) i_tmap 
                  WRITE ( *, 1100) nc, nr 
               ENDIF 
            ENDIF 
!                                                                       
!----- --- Type CON*                                                    
!                                                                       
         ELSEIF (icons.gt.0) then 
            READ (line (5:icons + 3), *, err = 998) ib, nc, nr 
            lbank = (ib.eq.bid (ibank) ) 
            IF (lbank) then 
               bank_line = line 
               READ (line (icons + 9:ll), *, err = 998) c_min, c_d 
               IF (dbg) then 
                  WRITE ( *, 1200) c_min, c_d 
                  WRITE ( *, 1100) nc, nr 
               ENDIF 
            ENDIF 
!                                                                       
!----- --- Type SLOG                                                    
!                                                                       
         ELSEIF (ilog.gt.0) then 
            READ (line (5:ilog + 3), *, err = 998) ib, nc, nr 
            lbank = (ib.eq.bid (ibank) ) 
            IF (lbank) then 
               bank_line = line 
               READ (line (ilog + 9:ll), *, err = 998) gtmin, gtmax,    &
               gtlog                                                    
               IF (dbg) then 
                  WRITE ( *, 1250) gtmin, gtmax, gtlog 
                  WRITE ( *, 1100) nc, nr 
               ENDIF 
            ENDIF 
!                                                                       
!----- --- Type unknown                                                 
!                                                                       
         ELSE 
            ier_num = - 48 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDIF 
      IF (.not.lbank) goto 20 
!                                                                       
!----- Now we can read the data                                         
!                                                                       
      IF (scal.ne.1.0) then 
         WRITE (output_io, 2000) ib, scal 
      ENDIF 
!                                                                       
      IF (index (bank_line, 'STD') .ne.0) then 
         READ (ifil, '(10(i2,f6.0))') (idummy, yval (i) , i = 1, nc) 
         DO i = 1, nc 
         yval (i) = scal * yval (i) 
         dyval (i) = sqrt (yval (i) ) 
         ENDDO 
         IF (dbg) WRITE ( * , 1300) 'STD' 
!                                                                       
      ELSEIF (index (bank_line, 'ESD') .ne.0) then 
         read (ifil, '(10f8.0)') (yval (i) , dyval (i) , i = 1, nc)
         DO i = 1, nc 
         yval (i) = scal * yval (i) 
         dyval (i) = scal * dyval (i) 
         ENDDO 
         IF (dbg) WRITE ( * , 1300) 'ESD' 
!                                                                       
      ELSEIF (index (bank_line, 'FXYE') .ne.0) then 
         DO i = 1, nc 
         READ (ifil, * ) xval (i), yval (i), dyval (i) 
         yval (i) = scal * yval (i) 
         dyval (i) = scal * dyval (i) 
         ENDDO 
         IF (dbg) WRITE ( * , 1300) 'FXYE' 
!                                                                       
      ELSE 
         IF (dbg) WRITE ( * , 1300) '(default) STD' 
         READ (ifil, '(10(i2,f6.0))') (idummy, yval (i) , i = 1, nc) 
         DO i = 1, nc 
         yval (i) = scal * yval (i) 
         dyval (i) = sqrt (yval (i) ) 
         ENDDO 
      ENDIF 
      ndat = nc 
!                                                                       
!----- Now we generate the TOF entry - ICONS                            
!                                                                       
      IF (icons.gt.0) then 
         DO i = 1, nc 
         xval (i) = c_min + (i - 1) * c_d 
         yval (i) = yval (i) / c_d 
         dyval (i) = dyval (i) / c_d 
         ENDDO 
         tof_offset = 0.0 
!                                                                       
!----- SLOG                                                             
!                                                                       
      ELSEIF (ilog.gt.0) then 
         DO i = 1, nc - 1 
         xdelta = xval (i + 1) - xval (i) 
         xval (i) = xval (i) + 0.5 * xdelta 
         yval (i) = yval (i) / 10.0 / xdelta 
         dyval (i) = dyval (i) / 10.0 / xdelta 
         ENDDO 
         xval (nc) = gtmax 
         xdelta = xval (nc) - xval (nc - 1) 
         yval (nc) = yval (nc) / 10.0 / xdelta 
         dyval (nc) = dyval (nc) / 10.0 / xdelta 
!                                                                       
!----- TIME_MAP                                                         
!                                                                       
      ELSEIF (itmap.gt.0) then 
!                                                                       
!----- we read the first time_map anyway (due to a bug in FSTBUSBIN)    
!----- then try finding the matching one                                
!----- (HIPPO will have mult. TIME_MAPs)                                
!                                                                       
         REWIND (ifil) 
   50    CONTINUE 
         READ (ifil, '(a)', end = 60, err = 998) line 
         ll = len_str (line) 
         IF (line (1:8) .eq.'TIME_MAP') then 
            lt = index (line (9:ll) , 'TIME_MAP') + 8 
            READ (line (9:lt), *, err = 998) it_no, it_nval, it_nrec 
            READ (line (lt + 9:ll), *, err = 998) clckwdt 
            maxmapu = (it_nval - 1) / 3 
            READ (ifil, '(10i8)', err = 998, end = 998) ( (tmap (j, i) ,&
            j = 1, 3) , i = 1, maxmapu) , ntmax                         
            ltmap = .true. 
         ENDIF 
         IF (.not.ltmap.or.it_no.ne.i_tmap) goto 50 
   60    CONTINUE 
!                                                                       
         IF (ltmap) then 
            IF (dbg) then 
               WRITE ( *, 1400) it_no, it_nval, it_nrec 
               WRITE ( *, 1500) clckwdt 
            ENDIF 
            IF (it_no.ne.i_tmap.and.dbg) WRITE ( *, 1600) 
         ELSE 
            ier_num = - 49 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
!----- - Now generate the tof entry from the tmap array                 
!                                                                       
         maxmapu = (it_nval - 1) / 3 
!                                                                       
         iq = 0 
         tof_offset = REAL(tmap (3, 1) ) / 20.0 
         DO i = 1, maxmapu - 1 
         itmap3i = tmap (3, i) 
         tmap3i = tmap (3, i) 
         itmap2i = tmap (2, i) 
         itmap1i = tmap (1, i) 
         itmap1iplus1 = tmap (1, i + 1) 
         DO j = itmap1i, itmap1iplus1 - 1 
         iq = iq + 1 
         xval (iq) = itmap2i + itmap3i * (j - itmap1i) + tmap3i / 2.0 
         yval (iq) = yval (iq) / itmap3i 
         dyval (iq) = dyval (iq) / itmap3i 
         ENDDO 
         ENDDO 
                                                                        
         itmap3i = tmap (3, maxmapu) 
         tmap3i = tmap (3, maxmapu) 
         itmap2i = tmap (2, maxmapu) 
         itmap1i = tmap (1, maxmapu) 
         imax = (ntmax - itmap2i) / itmap3i + itmap1i 
         DO j = itmap1i, imax 
         iq = iq + 1 
         xval (iq) = itmap2i + itmap3i * (j - itmap1i) + tmap3i / 2.0 
         yval (iq) = yval (iq) / itmap3i 
         dyval (iq) = dyval (iq) / itmap3i 
         ENDDO 
         onetick = clckwdt / 1000. 
         DO iq = 1, imax 
         xval (iq) = xval (iq) * onetick 
         ENDDO 
      ENDIF 
      ndat = nc - 1 
      RETURN 
!                                                                       
  510 CONTINUE 
      ier_num = - 50 
      ier_typ = ER_APPL 
      RETURN 
!                                                                       
  998 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
 1000 FORMAT    (' debug > TOF binning found: TIME_MAP number ',i3) 
 1100 FORMAT    (' debug > Bank has ',i6,' channels and ',i6,' records') 
 1200 FORMAT    (' debug > TOF binning CONST found: TOF_MIN: ',f9.2,    &
     &                  ', TOF_D: ',f9.2)                               
 1250 FORMAT    (' debug > TOF binning SLOG found:  TOF_MIN: ',f9.2,    &
     &                  ' TOF_MIN: ',f9.2,' DT/T: ',f9.6)               
 1300 FORMAT    (' debug > Data type ',a,' found') 
 1400 FORMAT    (' debug > TIME_MAP found: no: ',i3,', nval: ',i6,      &
     &                  ', nrec: ',i6)                                  
 1500 FORMAT    (' debug > TIME_MAP found: Clockwidth: ',f10.4) 
 1600 FORMAT    (' debug > No matching TIME_MAP found, taking first !') 
 2000 FORMAT    (' Rescaling bank ',i3,' by ',f8.1) 
 5000 FORMAT    (' Monitor count found: ',f10.1) 
!                                                                       
      END SUBROUTINE read_gsas_data                 
!*****7**************************************************************** 
      SUBROUTINE read_iparm (inam, nbank, difc, difa, tzero, twotheta,  &
      vtype, vrange, vp)                                                
!+                                                                      
!     Read GSAS instrument parameter file                               
!-                                                                      
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
USE lib_length
USE support_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER nbank 
!                                                                       
      CHARACTER ( * ) inam 
      CHARACTER(80) line 
      CHARACTER(6) search 
      REAL(8) difc (nbank) 
      REAL(8) difa (nbank) 
      REAL(8) tzero (nbank) 
      REAL vrange (nbank, 2), vp (nbank, 12) 
      REAL twotheta (nbank) 
      REAL dummy 
      INTEGER vtype (nbank) 
      INTEGER ibank, i, ll 
!                                                                       
!                                                                       
      CALL oeffne (12, inam, 'old') 
      IF (ier_num.ne.0) return 
!                                                                       
!------ Read instrument parameter file information                      
!                                                                       
   20 CONTINUE 
      READ (12, '(a)', end = 40) line 
      ll = len_str (line) 
      IF (ll.eq.0) goto 20 
      READ (line (1:ll) , '(4x,i2,a6)', err = 20) ibank, search 
!                                                                       
!------ Diffractometer constants                                        
!                                                                       
      IF (search.eq.' ICONS') then 
         IF (ibank.le.0.or.ibank.gt.nbank) goto 997 
         READ (line (13:ll), *, err = 998, end = 998) difc (ibank),     &
         difa (ibank), tzero (ibank)                                    
!                                                                       
!------ Bank parameters                                                 
!                                                                       
      ELSEIF (search.eq.'BNKPAR') then 
         IF (ibank.le.0.or.ibank.gt.nbank) goto 997 
         READ (line (13:ll), *, err = 998) dummy, twotheta (ibank) 
!                                                                       
!------ Incident spectrum information                                   
!                                                                       
      ELSEIF (search.eq.'I ITYP') then 
         IF (ibank.le.0.or.ibank.gt.nbank) goto 997 
         READ (line (13:ll), *, err = 998) vtype (ibank), vrange (ibank,&
         1), vrange (ibank, 2), dummy                                   
      ELSEIF (search.eq.'ICOFF1') then 
         IF (ibank.le.0.or.ibank.gt.nbank) goto 997 
         READ (line (13:ll), *, err = 998) (vp (ibank, i), i = 1, 4) 
      ELSEIF (search.eq.'ICOFF2') then 
         IF (ibank.le.0.or.ibank.gt.nbank) goto 997 
         READ (line (13:ll), *, err = 998) (vp (ibank, i), i = 5, 8) 
      ELSEIF (search.eq.'ICOFF3') then 
         IF (ibank.le.0.or.ibank.gt.nbank) goto 997 
         READ (line (13:ll), *, err = 998) (vp (ibank, i), i = 9, 12) 
!                                                                       
      ENDIF 
      GOTO 20 
!                                                                       
   40 CONTINUE 
      CLOSE (12) 
!                                                                       
      IF (dbg) then 
         WRITE (output_io, 1000) 
         DO i = 1, nbank 
         IF (twotheta (i) .ne. - 9999.) then 
            WRITE (output_io, 1100) i, twotheta (i), difc (i), difa (i),&
            tzero (i)                                                   
         ENDIF 
         ENDDO 
         WRITE (output_io, 1200) 
      ENDIF 
!                                                                       
      RETURN 
!                                                                       
  997 CONTINUE 
      ier_num = - 46 
      ier_typ = ER_APPL 
      CLOSE (12) 
      RETURN 
!                                                                       
  998 CONTINUE 
      ier_num = - 47 
      ier_typ = ER_APPL 
      CLOSE (12) 
!                                                                       
 1000 FORMAT (' debug > ',49('-'),/,                                    &
     &        ' debug > Bank',4x,'TwoTheta',5x,'DifC',9x,'DifA',        &
     &        6x,'Tzero')                                               
 1100 FORMAT (' debug > ',i3,5x,f7.2,3x,f9.2,3x,f8.2,3x,f8.2) 
 1200 FORMAT (' debug > ',49('-')) 
      END SUBROUTINE read_iparm                     
!*****7**************************************************************** 
      INTEGER function gsas_no_banks (ifil, bid, nbank) 
!+                                                                      
!     Determined number of banks in GSAS file                           
!-                                                                      
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
USE lib_length
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ifil, nbank 
      INTEGER bid (nbank) 
!                                                                       
      CHARACTER(80) line 
      INTEGER ibank, ib 
!                                                                       
!                                                                       
      ibank = 0 
      gsas_no_banks = 0
   10 CONTINUE 
      READ (ifil, '(a)', end = 20) line 
      IF (line (1:5) .eq.'BANK ') then 
         READ (line (6:len_str (line) ), *, err = 998) ib 
         IF (ib.ne.0) then 
            ibank = ibank + 1 
            bid (ibank) = ib 
         ENDIF 
      ENDIF 
      GOTO 10 
   20 CONTINUE 
!                                                                       
      gsas_no_banks = ibank 
      IF (dbg) WRITE ( *, 1000) ibank 
!                                                                       
      REWIND (ifil) 
      RETURN 
!                                                                       
  998 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
 1000 FORMAT(' debug > Number of GSAS banks found: ',i3) 
      END FUNCTION gsas_no_banks                    
!*****7**************************************************************** 
!
SUBROUTINE read_xy (ifil, ianz, werte, maxw) 
!+                                                                      
!     Load x,y or x,y,dx,dy files                                       
!-                                                                      
USE blanks_mod
USE errlist_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
USE precision_mod
!
USE count_col_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: mm = 50
!                                                                       
INTEGER, INTENT(IN) :: ifil 
INTEGER, INTENT(IN) :: ianz 
INTEGER, INTENT(IN) :: MAXW 
REAL(KIND=PREC_DP), INTENT(IN) :: werte (maxw) 
!                                                                       
CHARACTER(LEN=PREC_STRING) :: line 
REAL    :: values (mm) 
INTEGER :: nr, nval, iwex, iwey, iwdx, iwdy, iski 
INTEGER :: i, maxpp 
INTEGER :: ios      ! I/O status
INTEGER :: length   ! Character string length
!                                                                       
!------ Get correct column numbers                                      
!                                                                       
IF (ianz.eq.2) then 
   iwex = nint (werte (1) ) 
   iwey = nint (werte (2) ) 
   iwdx = 0 
   iwdy = 0 
   iski = 0 
ELSEIF (ianz.eq.3) then 
   iwex = nint (werte (1) ) 
   iwey = nint (werte (2) ) 
   iwdx = 0 
   iwdy = nint (werte (3) ) 
   iski = 0 
ELSEIF (ianz.eq.4) then 
   iwex = nint (werte (1) ) 
   iwey = nint (werte (2) ) 
   iwdx = nint (werte (3) ) 
   iwdy = nint (werte (4) ) 
   iski = 0 
ELSEIF (ianz.eq.5) then 
   iwex = nint (werte (1) ) 
   iwey = nint (werte (2) ) 
   iwdx = nint (werte (3) ) 
   iwdy = nint (werte (4) ) 
   iski = nint (werte (5) ) 
ELSE 
   iski = 0 
   iwex = 1 
   iwey = 2 
   iwdx = 0 
   iwdy = 0 
ENDIF 
!                                                                       
!------ read data                                                       
!                                                                       
nr = 1 
maxpp = maxarray - offxy (iz - 1) 
!                                                                       
ios = 0
DO i = 1, iski 
!  READ (ifil, 9999, end = 20) line 
   READ(ifil,'(a)', IOSTAT=ios ) line
   IF(IS_IOSTAT_END(ios)) RETURN
ENDDO 
!                                                                       
!  10 CONTINUE 
loop_col: DO
   READ(ifil,'(a)', IOSTAT=ios ) line
   IF(IS_IOSTAT_END(ios)) RETURN
   length = LEN_TRIM(line)
   CALL rem_leading_bl (line, length)
!  IF (line (1:1) .eq.'#'.and.iski.eq.0) goto 10 
   IF (line (1:1) .eq.'#'.and.iski.eq.0) CYCLE loop_col
   CALL count_col (line, nval) 
!  IF (nval.eq.0) goto 10 
   IF (nval >  0) EXIT loop_col
enddo loop_col
!
IF (nval.gt.mm) nval = mm 
BACKSPACE (ifil) 
!                                                                       
IF (ianz.lt.2) then 
   IF (nval.eq.3) iwdy = 3 
   IF (nval.ge.4) then 
      iwdx = 3 
      iwdy = 4 
   ENDIF 
ENDIF 
!                                                                       
IF (iski.gt.0) WRITE (output_io, 1010) iski 
WRITE (output_io, 1000) iwex, iwey, iwdx, iwdy 
IF (iwex.le.0 .OR. iwex.gt.nval .OR. iwey.le.0 .OR. iwey.gt.nval .OR. &
    iwdx.gt.nval .OR. iwdy.gt.nval) THEN                                    
   ier_num = - 59 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
ios = 0
loop_read: DO
   READ(ifil,'(a)', IOSTAT=ios ) line
   IF(IS_IOSTAT_END(ios)) THEN
      ios = 0
      EXIT loop_read
   ENDIF
   length = LEN_TRIM(line)
   CALL rem_leading_bl (line, length)
!  IF(line(1:1)=='#') GOTO 15
   IF(line(1:1)=='#') CYCLE loop_read
   READ (line, *, IOSTAT=ios, ERR=50 ) (values (i), i = 1, nval) 
   IF(IS_IOSTAT_END(ios)) THEN
      CYCLE loop_read
   ELSEIF(ios/=0) THEN
      EXIT loop_read
   ENDIF
   x (offxy (iz - 1) + nr) = values (iwex) 
   y (offxy (iz - 1) + nr) = values (iwey) 
   dx (offxy (iz - 1) + nr) = 0.0 
   dy (offxy (iz - 1) + nr) = 0.0 
   IF (iwdx.gt.0) dx (offxy (iz - 1) + nr) = values (iwdx) 
   IF (iwdy.gt.0) dy (offxy (iz - 1) + nr) = values (iwdy) 
   nr = nr + 1 
   IF (nr.gt.maxpp) then 
      ier_num = - 6 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
ENDDO loop_read
!
!20 CONTINUE 
IF(ios==0) THEN  
!
   lenc(iz) = nr - 1 
   offxy (iz) = offxy (iz - 1) + lenc(iz) 
   offz (iz) = offz (iz - 1) 
   iz = iz + 1 
!                                                                       
   CALL show_data (iz - 1) 
ENDIF
!
50 CONTINUE
!
IF(ios /= 0) THEN
   ier_num = -3
   ier_typ = ER_IO
ENDIF
!                                                                       
 1000 FORMAT     (1x,'Column assignement (x y dx dy): ',4i5) 
 1010 FORMAT     (1x,'Header lines skipped: ',i5) 
 9999 FORMAT     (a) 
!                                                                       
      END SUBROUTINE read_xy                        
!*****7**************************************************************** 
SUBROUTINE do_read_csv(MAXW, ianz, cpara, lpara, NOPTIONAL, opara, &
                       lopara, owerte, ncalc, ifil, lecho )
!
USE blanks_mod
USE charact_mod
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: O_SKIP      = 1
INTEGER, PARAMETER :: O_COLX      = 2
INTEGER, PARAMETER :: O_COLY      = 3
INTEGER, PARAMETER :: O_COLDX     = 4
INTEGER, PARAMETER :: O_COLDY     = 5
!INTEGER, PARAMETER :: O_LAYER     = 5
INTEGER, PARAMETER :: O_SEPARATOR = 7
INTEGER                                   , INTENT(IN   ) :: MAXW
INTEGER                                   , INTENT(INOUT) :: ianz
CHARACTER (LEN=PREC_STRING), DIMENSION(MAXW)     , INTENT(INOUT) :: cpara
INTEGER             , DIMENSION(MAXW)     , INTENT(INOUT) :: lpara
INTEGER                                   , INTENT(IN   ) :: NOPTIONAL
CHARACTER (LEN=PREC_STRING), DIMENSION(NOPTIONAL), INTENT(INOUT) :: opara
INTEGER             , DIMENSION(NOPTIONAL), INTENT(INOUT) :: lopara
REAL(KIND=PREC_DP)  , DIMENSION(NOPTIONAL), INTENT(INOUT) :: owerte
INTEGER                                   , INTENT(IN   ) :: ncalc
INTEGER                                   , INTENT(IN   ) :: ifil
LOGICAL                                   , INTENT(IN   ) :: lecho
!
CHARACTER(LEN=PREC_STRING) :: line
INTEGER             :: ios
INTEGER :: isep, iskip, icolx, icoly, icoldx, icoldy, ncol
INTEGER :: length
INTEGER :: nsep, nsec
INTEGER :: i, j
INTEGER             :: nr
INTEGER             :: maxpp
!INTEGER, DIMENSION(2,4) :: sections
INTEGER, DIMENSION(:,:), ALLOCATABLE :: limits
!
! Interpret the optional parameters into local variables
!
isep = 0
IF(opara(O_SEPARATOR)==';' .OR. opara(O_SEPARATOR)=='semicolon' ) THEN
   isep = IACHAR(';')
ELSEIF(opara(O_SEPARATOR)==':' .OR. opara(O_SEPARATOR)=='colon' ) THEN
   isep = IACHAR(':')
ELSEIF(opara(O_SEPARATOR)=='comma' ) THEN
   isep = IACHAR(',')
ELSEIF(opara(O_SEPARATOR)(1:3)=='tab' ) THEN
   isep = 9
ELSEIF(opara(O_SEPARATOR)(1:5)=='blank' ) THEN
   isep = blank1
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = 'Unknown separator'
ENDIF
!
iskip  = NINT(owerte(O_SKIP))
icolx  = NINT(owerte(O_COLX))
icoly  = NINT(owerte(O_COLY))
icoldx = NINT(owerte(O_COLDX))
icoldy = NINT(owerte(O_COLDY))
ncol   = MAX(icolx, icoly, icoldx, icoldy)
!
!write(*,*) ' IANZ ', ianz
!do i=1,ianz
!   write(*,*) ' PARA ', i, cpara(i)(1:lpara(i))
!enddo
!write(*,*) ' IANZ ', NOPTIONAL
!do i=1,NOPTIONAL
!   write(*,*) ' PARA ', i, opara(i)(1:lopara(i))
!enddo
!
!  Skip header lines
!
skip: DO i=1, iskip
   READ(ifil,'(a)', IOSTAT=ios) line
   IF(ios/=0) RETURN
ENDDO skip
!
! Main body read lines, determin the number of sections and read 
! the intended sections into x, y, dx, dy
! If any of the requested columns is empty stop reading
!
nr = 1
maxpp = maxarray - offxy (iz - 1) 
body: DO
   READ(ifil,'(a)', IOSTAT=ios) line
   length = LEN_TRIM(line)
   IF(ios/=0) EXIT body             ! Error, finish read
   IF(isep==blank1) THEN
      CALL rem_leading_bl(line,length)
      CALL rem_dbl_bl(line,length)
   ENDIF
   nsep = 0
   DO i=1,length                    ! Determin number of separators
      if(IACHAR(line(i:i))==isep) nsep = nsep + 1
   ENDDO
   nsec = nsep
!                                   ! last char /= sep => one more section
   IF(IACHAR(line(length:length))/=isep) nsec = nsec + 1
   IF(nsec<ncol) EXIT body          ! Less sections than required columns
   ALLOCATE(limits(2,nsec))
   limits(:,:)    = 0
   limits(1,1) = 1
   limits(1,2) = 1
   limits(1,nsec) = 1
   limits(2,nsec) = length
   j = 1                            ! Determine the limits of each section
   DO i=1, length
      IF(IACHAR(line(i:i))==isep) THEN    ! Next section boundary
         limits(2,j) = i-1
         limits(1,j+1) = MIN(i+1, length)
         j = j+1
      ENDIF
   ENDDO
!
!  Check if we have something like ...;;... for a required column If so quit reading
   IF(limits(2,icolx)-limits(1,icolx)<0) EXIT body   ! Width is zero 
   IF(limits(2,icoly)-limits(1,icoly)<0) EXIT body   ! Width is zero 
   IF(icoldx>0) THEN
      IF( limits(2,icoldx)-limits(1,icoldx)<0) EXIT body   ! Width is zero 
   ENDIF
   IF(icoldy>0) THEN
      IF( limits(2,icoldy)-limits(1,icoldy)<0) EXIT body   ! Width is zero 
   ENDIF
!
!  Read data x, y, optionally dx, dy
   READ(line(limits(1,icolx):limits(2,icolx)),*,IOSTAT=ios) x (offxy (iz - 1) + nr)
   IF(ios/=0) EXIT body
   READ(line(limits(1,icoly):limits(2,icoly)),*,IOSTAT=ios) y (offxy (iz - 1) + nr)
   IF(ios/=0) EXIT body
   IF(icoldx>0) THEN
      READ(line(limits(1,icoldx):limits(2,icoldx)),*,IOSTAT=ios) dx (offxy (iz - 1) + nr)
      IF(ios/=0) EXIT body
   ENDIF
   IF(icoldy>0) THEN
      READ(line(limits(1,icoldy):limits(2,icoldy)),*,IOSTAT=ios) dy (offxy (iz - 1) + nr)
      IF(ios/=0) EXIT body
   ENDIF
   nr = nr + 1                            ! Increment number of data points
   IF (nr.gt.maxpp) THEN                  ! Too many data points ERROR
      ier_num = - 6 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
   DEALLOCATE(limits)
!
ENDDO body
!
! Record data set length and show data
!
lenc(iz)   = nr - 1 
offxy(iz) = offxy(iz - 1) + lenc(iz) 
offz(iz)  = offz(iz - 1) 
iz = iz + 1 
!                                                                       
IF(lecho) CALL show_data (iz - 1) 
!                                                                       
END SUBROUTINE do_read_csv
!*****7**************************************************************** 
      SUBROUTINE read_csv (ifil)
!+                                                                      
!     Load CSV Files generated from PANAlytical files                   
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!
      INTEGER, INTENT(in) :: ifil
!
      CHARACTER(LEN=PREC_STRING) :: line
      CHARACTER(LEN=1   ) :: sep
      INTEGER             :: ios
      INTEGER             :: isep
      INTEGER             :: length
      INTEGER             :: nr
      INTEGER             :: maxpp
!
      header: DO
         READ(ifil, '(a)', IOSTAT=ios) line
         IF(ios/=0) RETURN
         IF(line(1:5) == 'Angle') THEN
            sep = line(6:6)
            EXIT header
         ENDIF
      ENDDO header
!
      nr = 1
      maxpp = maxarray - offxy (iz - 1) 
      body: DO
         READ(ifil, '(a)', IOSTAT=ios) line
         IF(ios/=0) EXIT body
         isep = INDEX(line,sep)
         length = LEN_TRIM(line)
         READ(line(1:isep-1),*,IOSTAT=ios)      x (offxy (iz - 1) + nr)
         READ(line(isep+1:length),*,IOSTAT=ios) y (offxy (iz - 1) + nr)
         dx(offxy (iz - 1) + nr) = 0.0
         dy(offxy (iz - 1) + nr) = 0.0
         nr = nr + 1
         IF (nr.gt.maxpp) then 
            ier_num = - 6 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDDO body
!
      lenc (iz) = nr - 1 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
      offz (iz) = offz (iz - 1) 
      iz = iz + 1 
!                                                                       
      CALL show_data (iz - 1) 
!                                                                       
      END SUBROUTINE read_csv
!*****7**************************************************************** 
      SUBROUTINE read_mca_single (ifil, ianz, werte, maxw) 
!+                                                                      
!     Load single MCA spectrum                                          
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      REAL(KIND=PREC_DP) :: werte (maxw)
      REAL :: xold, xnew, ynew 
      INTEGER nr, ifil, ianz, iscan, nscan, maxpp 
!                                                                       
!------ check parameters                                                
!                                                                       
      IF (ianz.ne.1) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
!------ Move forward to correct scan                                    
!                                                                       
      iscan = nint (werte (1) ) 
      nscan = 1 
      xold = - 9999. 
!                                                                       
    5 CONTINUE 
      IF (nscan.eq.iscan) then 
         BACKSPACE (ifil) 
         GOTO 7 
      ENDIF 
      READ (ifil, *, end = 90) xnew, ynew 
      IF (xnew.lt.xold) nscan = nscan + 1 
      xold = xnew 
      GOTO 5 
    7 CONTINUE 
!                                                                       
!------ read data                                                       
!                                                                       
      nr = 1 
      xold = - 9999. 
!                                                                       
      maxpp = maxarray - offxy (iz - 1) 
   10 CONTINUE 
      READ (ifil, *, end = 20) xnew, ynew 
      IF (xnew.lt.xold) goto 20 
      x (offxy (iz - 1) + nr) = xnew 
      y (offxy (iz - 1) + nr) = ynew 
      dx (offxy (iz - 1) + nr) = 0.0 
      dy (offxy (iz - 1) + nr) = 0.0 
      xold = xnew 
      nr = nr + 1 
      IF (nr.gt.maxpp) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      GOTO 10 
   20 CONTINUE 
      lenc (iz) = nr - 1 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
      offz (iz) = offz (iz - 1) 
      fform (iz) = 'XY' 
      iz = iz + 1 
!                                                                       
      WRITE (titel (iwin, iframe, 2), 1000) iscan 
      achse (iwin, iframe, 1) = 'Channel' 
      achse (iwin, iframe, 2) = 'Counts' 
!                                                                       
      CALL show_data (iz - 1) 
      RETURN 
!                                                                       
   90 CONTINUE 
      ier_num = - 41 
      ier_typ = ER_APPL 
!                                                                       
 1000 FORMAT     ('MCA spectrum #',i5) 
!                                                                       
      END SUBROUTINE read_mca_single                
!*****7**************************************************************** 
      SUBROUTINE read_xx (ifil) 
!+                                                                      
!     Load files point no versus value                                  
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER ifil, nr 
      INTEGER maxpp 
!                                                                       
      nr = 1 
      maxpp = maxarray - offxy (iz - 1) 
   10 CONTINUE 
      READ (ifil, *, end = 20) y (offxy (iz - 1) + nr) 
      x (offxy (iz - 1) + nr) = REAL(nr) 
      dx (offxy (iz - 1) + nr) = 0.0 
      dy (offxy (iz - 1) + nr) = 0.0 
      nr = nr + 1 
      IF (nr.gt.maxpp) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      GOTO 10 
   20 CONTINUE 
      lenc (iz) = nr - 1 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
      offz (iz) = offz (iz - 1) 
      fform (iz) = 'XY' 
      iz = iz + 1 
!                                                                       
      CALL show_data (iz - 1) 
!                                                                       
      END SUBROUTINE read_xx                        
!*****7**************************************************************** 
      SUBROUTINE read_crystal (ifil) 
!+                                                                      
!     Reads special crystal file format exported from DISCUS.           
!     Line per atom  : x,y,z, mtyp,mcol,msiz                            
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL rmsiz, zdummy 
      INTEGER ifil, imcol, imtyp, nr, iii, maxpp 
!                                                                       
      nr = 1 
      maxpp = maxarray - offxy (iz - 1) 
   10 CONTINUE 
      iii = offxy (iz - 1) + nr 
      READ (ifil, *, end = 20) x (iii), y (iii), zdummy, imtyp, imcol,  &
      rmsiz                                                             
!                                                                       
!-------  Put color / marker and size information in DX, DY             
!                                                                       
      dx (iii) = REAL(imcol * 1000 + imtyp) 
      dy (iii) = rmsiz 
!                                                                       
      nr = nr + 1 
      IF (nr.gt.maxpp) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      GOTO 10 
   20 CONTINUE 
      lenc (iz) = nr - 1 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
      offz (iz) = offz (iz - 1) 
      iz = iz + 1 
!                                                                       
      CALL show_data (iz - 1) 
!                                                                       
      END SUBROUTINE read_crystal                   
!*****7**************************************************************** 
      SUBROUTINE read_marker (ifil) 
!+                                                                      
!     Read marker file (only x-values)                                  
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: line 
      INTEGER ifil, nr, maxpp 
!                                                                       
!------ read data                                                       
!                                                                       
      nr = 1 
      maxpp = maxarray - offxy (iz - 1) 
   10 CONTINUE 
      READ (ifil, 9999, end = 20) line 
      IF (line (1:1) .eq.'#') goto 10 
      BACKSPACE (ifil) 
      READ (ifil, *, end = 20) x (offxy (iz - 1) + nr) 
      y (offxy (iz - 1) + nr) = 0.0 
      dx (offxy (iz - 1) + nr) = 0.0 
      dy (offxy (iz - 1) + nr) = 0.0 
      nr = nr + 1 
      IF (nr.gt.maxpp) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      GOTO 10 
   20 CONTINUE 
      lenc (iz) = nr - 1 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
      offz (iz) = offz (iz - 1) 
      iz = iz + 1 
!                                                                       
!------ Some smart defaults                                             
!                                                                       
      ilinetyp (iwin, iframe, iz - 1) = 0 
      imarktyp (iwin, iframe, iz - 1) = 9 
      sizemark (iwin, iframe, iz - 1) = 1.0 
!                                                                       
      CALL show_data (iz - 1) 
!                                                                       
 9999 FORMAT    (a) 
      END SUBROUTINE read_marker                    
!*****7**************************************************************** 
      SUBROUTINE read_z (ifil, izaehl) 
!+                                                                      
!     Read powder scan files MAN1 program ..                            
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL dummy 
      INTEGER ifil, izaehl, nr, maxpp 
!                                                                       
!------ setup                                                           
!                                                                       
      achse (iwin, iframe, 1) = '2-theta' 
      achse (iwin, iframe, 2) = 'counts' 
!                                                                       
!------ read data                                                       
!                                                                       
      nr = 1 
      maxpp = maxarray - offxy (iz - 1) 
!                                                                       
      READ (ifil, *, end = 20) 
      READ (ifil, *, end = 20) 
      READ (ifil, *, end = 20) 
   10 CONTINUE 
      IF (izaehl.eq.1) then 
         READ (ifil, *, end = 20) x (offxy (iz - 1) + nr), y (offxy (iz &
         - 1) + nr)                                                     
      ELSE 
         READ (ifil, *, end = 20) x (offxy (iz - 1) + nr), dummy, y (   &
         offxy (iz - 1) + nr)                                           
      ENDIF 
      dx (offxy (iz - 1) + nr) = 0.0 
      dy (offxy (iz - 1) + nr) = 0.0 
      nr = nr + 1 
      IF (nr.gt.maxpp) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      GOTO 10 
   20 CONTINUE 
      lenc (iz) = nr - 1 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
      offz (iz) = offz (iz - 1) 
      iz = iz + 1 
!                                                                       
      CALL show_data (iz - 1) 
!                                                                       
      END SUBROUTINE read_z                         
!*****7**************************************************************** 
      SUBROUTINE get_extrema 
!+                                                                      
!     get max/min values of files                                       
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ik 
!                                                                       
      DO ik = 1, iz - 1 
      IF (lni (ik) ) then 
         CALL get_extrema_xy (x, ik, nx (ik), xmin, xmax) 
         CALL get_extrema_xy (y, ik, ny (ik), ymin, ymax) 
         CALL get_extrema_z (z, ik, nx (ik), ny (ik), zmin, zmax) 
      ELSE 
         CALL get_extrema_xy (x, ik, lenc (ik), xmin, xmax) 
         CALL get_extrema_xy (y, ik, lenc (ik), ymin, ymax) 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE get_extrema                    
!*****7**************************************************************** 
      SUBROUTINE get_col (zeile, field, nf, colm) 
!+                                                                      
!     This subroutine splits string at spaces                           
!-                                                                      
USE lib_length
      IMPLICIT none 
!                                                                       
      INTEGER colm 
!                                                                       
      CHARACTER ( * ) zeile, field (0:colm) 
      INTEGER nf, i, is, ia 
      LOGICAL ein 
!                                                                       
!                                                                       
      nf = 0 
      ia = 1
      ein = .false. 
!                                                                       
      is = 1 
      DO while (zeile (is:is) .ne.' ') 
      is = is + 1 
      ENDDO 
!                                                                       
      DO i = is, len_str (zeile) 
      IF (zeile (i:i) .ne.' ') then 
         IF (.not.ein) ia = i 
         ein = .true. 
      ELSEIF (zeile (i:i) .eq.' '.and.ein) then 
         ein = .false. 
         nf = nf + 1 
         field (nf) = zeile (ia:i - 1) 
      ENDIF 
      ENDDO 
!                                                                       
      IF (ein) then 
         nf = nf + 1 
         field (nf) = zeile (ia:i - 1) 
      ENDIF 
!                                                                       
      END SUBROUTINE get_col                        
!*****7**************************************************************** 
      SUBROUTINE read_simref (ifil) 
!+                                                                      
!     Load *.plt file from SIMREF                                       
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: line 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER ifil, j, npkt 
      INTEGER ianz, i, maxpp 
      REAL tth_anf, tth_del, tth_end, temp 
!                                                                       
!                                                                       
!------ space left for additional data set ?                            
!                                                                       
      IF (iz + 4.gt.maxkurvtot) then 
         ier_num = - 1 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      ianz = maxw 
!                                                                       
!------ read data                                                       
!                                                                       
      maxpp = (maxarray - offxy (iz - 1) ) / 4 
!                                                                       
      READ (ifil, 9999, end = 20) line 
      READ (ifil, 9999, end = 20) line 
      READ (line (4:200), * ) tth_anf, tth_del, tth_end, temp, npkt 
      IF (npkt.gt.maxpp) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!     allocate offsets, lengths,names                                   
!                                                                       
      DO i = 1, 4 
      lenc (iz) = npkt 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
      offz (iz) = offz (iz - 1) 
      iz = iz + 1 
      ENDDO 
      fname (iz - 3) = 'calc' 
      fname (iz - 2) = 'diff' 
      fname (iz - 1) = 'back' 
!                                                                       
      DO i = 1, npkt 
      READ (ifil, *, end = 20) (werte (j), j = 1, ianz) 
      x (offxy (iz - 5) + i) = werte (1) 
      x (offxy (iz - 4) + i) = werte (1) 
      x (offxy (iz - 3) + i) = werte (1) 
      x (offxy (iz - 2) + i) = werte (1) 
      y (offxy (iz - 5) + i) = werte (2) 
      y (offxy (iz - 4) + i) = werte (3) 
      y (offxy (iz - 3) + i) = werte (4) 
      y (offxy (iz - 2) + i) = werte (5) 
      dx (offxy (iz - 5) + i) = 0.0 
      dx (offxy (iz - 4) + i) = 0.0 
      dx (offxy (iz - 3) + i) = 0.0 
      dx (offxy (iz - 2) + i) = 0.0 
      dy (offxy (iz - 5) + i) = 0.0 
      dy (offxy (iz - 4) + i) = 0.0 
      dy (offxy (iz - 3) + i) = 0.0 
      dy (offxy (iz - 2) + i) = 0.0 
      ENDDO 
   20 CONTINUE 
!                                                                       
      DO j = 4, 1, - 1 
      CALL show_data (iz - j) 
      ENDDO 
!                                                                       
 9999 FORMAT     (a) 
!                                                                       
      END SUBROUTINE read_simref                    
!*****7**************************************************************** 
      REAL function chan2kev (channel, mca_par) 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER channel 
      REAL mca_par (3) 
!                                                                       
      chan2kev = mca_par (1) * channel**2 + mca_par (2) * channel +     &
      mca_par (3)                                                       
!                                                                       
      END FUNCTION chan2kev                         
